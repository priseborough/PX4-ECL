#!/usr/bin/env python3

from sympy import *
from code_gen import *
import numpy as np

# q: quaternion describing rotation from frame 1 to frame 2
# returns a rotation matrix derived form q which describes the same
# rotation
def quat2Rot(q):
    q0 = q[0]
    q1 = q[1]
    q2 = q[2]
    q3 = q[3]

    Rot = Matrix([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
                  [2*(q1*q2 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 - q0*q1)],
                   [2*(q1*q3-q0*q2), 2*(q2*q3 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])

    return Rot

def create_cov_matrix(i, j):
    if j >= i:
        # return Symbol("P(" + str(i) + "," + str(j) + ")", real=True)
        # legacy array format
        return Symbol("P[" + str(i) + "][" + str(j) + "]", real=True)
    else:
        return 0

def quat_mult(p,q):
    r = Matrix([p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
                p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
                p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
                p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0]])

    return r

def create_symmetric_cov_matrix(n):
    # define a symbolic covariance matrix
    P = Matrix(n,n,create_cov_matrix)

    for index in range(n):
        for j in range(n):
            if index > j:
                P[index,j] = P[j,index]

    return P

def generate_code():
    print('Starting code generation:')
    print('Creating symbolic variables ...')

    dt = symbols("dt", real=True)  # dt
    g = symbols("g", real=True) # gravity constant

    r_hor_vel = symbols("R_hor_vel", real=True) # horizontal velocity noise variance
    r_ver_vel = symbols("R_vert_vel", real=True) # vertical velocity noise variance
    r_hor_pos = symbols("R_hor_pos", real=True) # horizontal position noise variance

    # inputs, integrated gyro measurements
    # delta angle x y z
    d_ang_x, d_ang_y, d_ang_z = symbols("dax day daz", real=True)  # delta angle x
    d_ang = Matrix([d_ang_x, d_ang_y, d_ang_z])

    # inputs, integrated accelerometer measurements
    # delta velocity x y z
    d_v_x, d_v_y, d_v_z = symbols("dvx dvy dvz", real=True)
    d_v = Matrix([d_v_x, d_v_y,d_v_z])

    u = Matrix([d_ang, d_v])

    # input noise
    d_ang_x_var, d_ang_y_var, d_ang_z_var = symbols("daxVar dayVar dazVar", real=True)

    d_v_x_var, d_v_y_var, d_v_z_var = symbols("dvxVar dvyVar dvzVar", real=True)

    var_u = Matrix.diag(d_ang_x_var, d_ang_y_var, d_ang_z_var, d_v_x_var, d_v_y_var, d_v_z_var)

    # define state vector

    # attitude quaternion
    qw, qx, qy, qz = symbols("q0 q1 q2 q3", real=True)
    q = Matrix([qw,qx,qy,qz])
    R_to_earth = quat2Rot(q)
    R_to_body = R_to_earth.T

    # velocity in NED local frame (north, east, down)
    vx, vy, vz = symbols("vn ve vd", real=True)
    v = Matrix([vx,vy,vz])

    # position in NED local frame (north, east, down)
    px, py, pz = symbols("pn pe pd", real=True)
    p = Matrix([px,py,pz])

    # delta angle bias x y z
    d_ang_bx, d_ang_by, d_ang_bz = symbols("dax_b day_b daz_b", real=True)
    d_ang_b = Matrix([d_ang_bx, d_ang_by, d_ang_bz])
    d_ang_true = d_ang - d_ang_b

    # delta velocity bias x y z
    d_vel_bx, d_vel_by, d_vel_bz = symbols("dvx_b dvy_b dvz_b", real=True)
    d_vel_b = Matrix([d_vel_bx, d_vel_by, d_vel_bz])
    d_vel_true = d_v - d_vel_b

    # earth magnetic field vector x y z
    ix, iy, iz = symbols("magN magE magD", real=True)
    i = Matrix([ix,iy,iz])

    # earth magnetic field bias in body frame
    ibx, iby, ibz = symbols("ibx iby ibz", real=True)

    ib = Matrix([ibx,iby,ibz])

    # wind in local NE frame (north, east)
    wx, wy = symbols("vwn, vwe", real=True)
    w = Matrix([wx,wy])

    # state vector at arbitrary time t
    state = Matrix([q, v, p, d_ang_b, d_vel_b, i, ib, w])

    print('Defining state propagation ...')
    # kinematic processes driven by IMU 'control inputs'
    q_new = quat_mult(q, Matrix([1, 0.5 * d_ang_true[0],  0.5 * d_ang_true[1],  0.5 * d_ang_true[2]]))
    v_new = v + R_to_earth * d_vel_true + Matrix([0,0,g]) * dt
    p_new = p + v * dt

    # static processes
    d_ang_b_new = d_ang_b
    d_vel_b_new = d_vel_b
    i_new = i
    ib_new = ib
    w_new = w

    # predicted state vector at time t + dt
    state_new = Matrix([q_new, v_new, p_new, d_ang_b_new, d_vel_b_new, i_new, ib_new, w_new])

    print('Computing state propagation jacobian ...')
    A = state_new.jacobian(state)
    G = state_new.jacobian(u)

    P = create_symmetric_cov_matrix(24)

    print('Computing covariance propagation ...')
    P_new = A * P * A.T + G * var_u * G.T

    for index in range(24):
        for j in range(24):
            if index > j:
                P_new[index,j] = 0

    print('Simplifying covariance propagation ...')
    P_new_simple = cse(P_new, symbols("PS0:400"), optimizations='basic')

    print('Writing covariance propagation to file ...')
    cov_code_generator = CodeGenerator("./generated/covariance_generated.cpp")
    cov_code_generator.print_string("Equations for covariance matrix prediction, without process noise!")
    cov_code_generator.write_subexpressions(P_new_simple[0])
    cov_code_generator.write_matrix(Matrix(P_new_simple[1]), "nextP", True, "[", "]")

    cov_code_generator.close()

    print('Code generation finished!')

if __name__ == "__main__":
    generate_code()
