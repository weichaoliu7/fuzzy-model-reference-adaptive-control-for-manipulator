#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <time.h>
#include "sine.h"
#include "cosine.h"
#include "inv_matrix.h"

// reference: [1]Golea N, Golea A, Benmahammed K. Fuzzy model reference adaptive control[J]. IEEE Transactions on Fuzzy Systems, 2002, 10(4): 436-444.
// [2]佟绍成,非线性系统的自适应模糊控制

// global variables declaration
#define PI 3.14159
#define ARRAY_SIZE 30000  // sampling times
#define n_joint 2         // number of robot manipulator joints

static double Ts = 0.001;  // sampling period
static double t0 = 0.0;    // start time
static double t1 = 30.0;   // end time
static double l;           // length of robot manipulator joints
static double m[n_joint];  // mass of robot manipulator joints
static double kv[n_joint]; // coefficient of viscous friction
static double kc[n_joint]; // coefficient of coulomb friction
static double g = 9.41;   // gravitational acceleration

// calculate matrix multiplication
void matrix_Multi(double *C, double *A, double *B, int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            *(C + j * cols2 + k) = 0.0;
            for (int g = 0; g < cols1; g++){
                *(C + j * cols2 + k) += *(A + j * cols1 + g) * *(B + g * cols2 + k);
            }
        }
    }
}

// calculate the transpose of matrix
void matrix_transpose(int rows, int cols, double matrix[rows][cols], double result[cols][rows]){
    for (int j = 0; j < rows; j++){
        for (int k = 0; k < cols; k++){
            result[k][j] = matrix[j][k];
        }
    }
}

// symbolic function
double sign(double a) {
    if (a > 0) {
        return 1.0;
    } else if (a < 0) {
        return -1.0;
    } else {
        return 0.0;
    }
}

struct _archive{
    double q1_archive[ARRAY_SIZE];
    double dq1_archive[ARRAY_SIZE];
    double q2_archive[ARRAY_SIZE];
    double dq2_archive[ARRAY_SIZE];
    double error1_archive[ARRAY_SIZE];
    double error2_archive[ARRAY_SIZE];
    double error1_velocity_archive[ARRAY_SIZE];
    double error2_velocity_archive[ARRAY_SIZE];
    double torque1_archive[ARRAY_SIZE];
    double torque2_archive[ARRAY_SIZE];
} archive;

Data q1_desired, dq1_desired, qd1_1, qd1_2, dqd1_1, dqd1_2;
Data q2_desired, dq2_desired, qd2_1, qd2_2, dqd2_1, dqd2_2;

struct Amp{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
};

struct M0{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
};

struct B0{
    double qd1_1, qd1_2;
    double dqd1_1, dqd1_2;
    double qd2_1, qd2_2;
    double dqd2_1, dqd2_2;
};

void SystemInput(Data *q1_desired, Data *qd1_1, Data *qd1_2, Data *dq1_desired, Data *dqd1_1, Data *dqd1_2, Data *q2_desired, Data *qd2_1, Data *qd2_2, Data *dq2_desired, Data *dqd2_1, Data *dqd2_2, double Ts, double t0, double t1){

    struct Amp amp; // amplitude
    amp.qd1_1 = PI;
    amp.qd1_2 = PI * 0.1;
    amp.dqd1_1 = PI * 0.5;
    amp.dqd1_2 = PI * 0.1 * 2;
    amp.qd2_1 = PI * 0.5;
    amp.qd2_2 = PI * 0.1;
    amp.dqd2_1 = PI * 0.5;
    amp.dqd2_2 = PI * 0.1 * 3;

    struct M0 m0; // angular frequency
    m0.qd1_1 = 0.5;
    m0.qd1_2 = 2;
    m0.dqd1_1 = 0.5;
    m0.dqd1_2 = 2;
    m0.qd2_1 = 1;
    m0.qd2_2 = 3;
    m0.dqd2_1 = 1;
    m0.dqd2_2 = 3;

    struct B0 b0; // vertical shift
    b0.qd1_1 = 0;
    b0.qd1_2 = 0;
    b0.dqd1_1 = 0;
    b0.dqd1_2 = 0;
    b0.qd2_1 = 0;
    b0.qd2_2 = 0;
    b0.dqd2_1 = 0;
    b0.dqd2_2 = 0;

    sine(qd1_1, Ts, t0, t1, amp.qd1_1, m0.qd1_1, b0.qd1_1);
    sine(qd1_2, Ts, t0, t1, amp.qd1_2, m0.qd1_2, b0.qd1_2);
    cosine(dqd1_1, Ts, t0, t1, amp.dqd1_1, m0.dqd1_1, b0.dqd1_1);
    cosine(dqd1_2, Ts, t0, t1, amp.dqd1_2, m0.dqd1_2, b0.dqd1_2);
    sine(qd2_1, Ts, t0, t1, amp.qd2_1, m0.qd2_1, b0.qd2_1);
    sine(qd2_2, Ts, t0, t1, amp.qd2_2, m0.qd2_2, b0.qd2_2);
    cosine(dqd2_1, Ts, t0, t1, amp.dqd2_1, m0.dqd2_1, b0.dqd2_1);
    cosine(dqd2_2, Ts, t0, t1, amp.dqd2_2, m0.dqd2_2, b0.dqd2_2);

    q1_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    q2_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    dq1_desired->y = malloc(sizeof(double) * ARRAY_SIZE);
    dq2_desired->y = malloc(sizeof(double) * ARRAY_SIZE);

    for (int j = 0; j < ARRAY_SIZE; j++){
        q1_desired->y[j] = qd1_1->y[j] + qd1_2->y[j];    // desired angular position of joint 1
        q2_desired->y[j] = qd2_1->y[j] + qd2_2->y[j];    // desired angular position of joint 2
        dq1_desired->y[j] = dqd1_1->y[j] + dqd1_2->y[j]; // desired angular velocity of joint 1
        dq2_desired->y[j] = dqd2_1->y[j] + dqd2_2->y[j]; // desired angular velocity of joint 2
    }
}

struct _system_state{
    double q[n_joint];   // actual angular displacement
    double dq[n_joint];  // actual angular velocity
    double ddq[n_joint]; // actual angular acceleration
} system_state;

struct _torque{
    double torque[n_joint];       // control input torque
    double torque_fuzzy[n_joint]; // output of fuzzy controller
    double torque_si[n_joint];    // additional control term used to overcome uncertainties
} torque;

struct _dynamics{
    double D[n_joint][n_joint]; // inertia matrix for manipulator dynamics equation
    double C[n_joint][n_joint]; // Coriolis and centrifugal matrix for manipulator dynamics equation
    double G[n_joint];          // gravitational torque matrix for manipulator dynamics equation
    double Fv[n_joint];         // viscous friction torque
    double Fc[n_joint];         // coulomb friction torque
} dynamics;

// membership function of input vectors in fuzzy sets
struct _membership{
    double q1[3]; // membership function of q1
    double q2[3]; // membership function of q2
} membership;

struct _controller{
    double controller_u[10];
    double controller_out[n_joint];
    double x1[n_joint];                                   // state variable of the first reference model
    double x2[n_joint];                                   // state variable of the second reference model
    double r[n_joint];                                    // bounded reference input of reference model
    double x[n_joint * 2];                                // state variable
    double xi[3][n_joint];                                // normalized firing strength vector, defined in Eq. 6
    double z[n_joint * 2 + 1][n_joint];                   // defined in Eq. 5
    double error[n_joint];                                // angular displacement error
    double error_velocity[n_joint];                       // angular velocity error
    double error_state[n_joint][n_joint];                 // tracking error of reference model, defined in Eq. 12
    double bc[n_joint][n_joint];                          // [0;...;1]
    double gamma[n_joint][n_joint];                       // fuzzy controller parameter update gain
    double P[n_joint][n_joint][n_joint];                  // solution of reference model Lyapunov equation
    double theta_derivative[3][n_joint * 2 + 1][n_joint]; // derivative of integral term of PI update law, Eq. 9
    double theta[3][n_joint * 2 + 1][n_joint];            // integral term of PI update law
    double phi[3][n_joint * 2 + 1][n_joint];              // proportional term of PI update law, Eq. 8
    double THETA[3][n_joint * 2 + 1][n_joint];            // fuzzy controller parameter update PI law, Eq. 37
    double b[n_joint][n_joint][2];                        // upper and lower bound of bii,smooth unknown functions, 0 means lower, 1 means upper
    double Beta[n_joint][n_joint];                        // known function
    double omega_upper[n_joint];                          // upper bound of minimum approximation error
    double eta_upper[n_joint];                            // upper bound of external disturbance
    double d_upper[n_joint];                              // upper bound of uncertain term
    double delta[n_joint];                                // parameter, positive constant
    double beta[n_joint];                                 // parameter, positive constants
} controller;

void CONTROLLER_init(){
    system_state.q[0] = 0.1;
    system_state.dq[0] = 2.0;
    system_state.q[1] = 0.1;
    system_state.dq[1] = 2.5;
    torque.torque[0] = 0.0;
    torque.torque[1] = 0.0;
    controller.controller_u[0] = q1_desired.y[0];
    controller.controller_u[1] = dq1_desired.y[0];
    controller.controller_u[2] = q2_desired.y[0];
    controller.controller_u[3] = dq2_desired.y[0];
    controller.controller_u[4] = system_state.q[0];
    controller.controller_u[5] = system_state.dq[0];
    controller.controller_u[6] = system_state.q[1];
    controller.controller_u[7] = system_state.dq[1];
    controller.controller_u[8] = torque.torque[0];
    controller.controller_u[9] = torque.torque[1];
    controller.gamma[1][0] = 20; // fuzzy controller parameter integral term update gain
    controller.gamma[1][1] = 500;
    controller.gamma[0][0] = 100; // fuzzy controller parameter proportional term update gain
    controller.gamma[0][1] = 10;
    controller.b[0][0][0] = 1.7; // lower bound of bii,smooth unknown functions
    controller.b[1][1][0] = 2;
    controller.b[0][1][1] = 1.3; // lower bound of bij,i != j,smooth unknown functions
    controller.b[1][0][1] = 1.3;
    controller.omega_upper[0] = 0.9; // upper bound of minimum approximation error
    controller.omega_upper[1] = 0.9;
    controller.eta_upper[0] = 0.78; // upper bound of external disturbance
    controller.eta_upper[1] = 3.13;

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            controller.bc[k][j] = k;
        }
    }

    for (int j = 0; j < n_joint; j++){
        controller.P[0][0][j] = 9.6875; // solution of reference model Lyapunov equation
        controller.P[0][1][j] = 0.4687;
        controller.P[1][0][j] = 0.4687;
        controller.P[1][1][j] = 0.3711;
    }

    for (int j = 0; j < n_joint; j++){
        controller.d_upper[j] = controller.omega_upper[j] + controller.eta_upper[j] / controller.b[j][j][0]; // upper bound of uncertain term
    }

    for (int j = 0; j < n_joint; j++){
        controller.delta[j] = 0.1;
        controller.beta[j] = 0.1;
    }
    
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                controller.theta[k][g][j] = 0.0; // integral term of parameter update PI law
            }
        }
    }
}

double CONTROLLER_realize(int i){
    controller.controller_u[0] = q1_desired.y[i];    // desired angular position of joint 1
    controller.controller_u[1] = dq1_desired.y[i];   // desired angular velocity of joint 1
    controller.controller_u[2] = q2_desired.y[i];    // desired angular position of joint 2
    controller.controller_u[3] = dq2_desired.y[i];   // desired angular velocity of joint 2
    controller.controller_u[4] = system_state.q[0];  // actual angular position of joint 1
    controller.controller_u[5] = system_state.dq[0]; // actual angular velocity of joint 1
    controller.controller_u[6] = system_state.q[1];  // actual angular position of joint 2
    controller.controller_u[7] = system_state.dq[1]; // actual angular velocity of joint 2
    controller.controller_u[8] = torque.torque[0];   // control input torque of joint 1, Eq. 11
    controller.controller_u[9] = torque.torque[1];   // control input torque of joint 2, Eq. 11
    archive.q1_archive[i] = controller.controller_u[4];
    archive.dq1_archive[i] = controller.controller_u[5];
    archive.q2_archive[i] = controller.controller_u[6];
    archive.dq2_archive[i] = controller.controller_u[7];

    controller.error[0] = q1_desired.y[i] - system_state.q[0];            // angular position tracking error of link 1
    controller.error_velocity[0] = dq1_desired.y[i] - system_state.dq[0]; // angular velocity tracking error of link 1
    controller.error[1] = q2_desired.y[i] - system_state.q[1];            // angular position tracking error of link 2
    controller.error_velocity[1] = dq2_desired.y[i] - system_state.dq[1]; // angular velocity tracking error of link 2
    archive.error1_archive[i] = controller.error[0];
    archive.error1_velocity_archive[i] = controller.error_velocity[0];
    archive.error2_archive[i] = controller.error[1];
    archive.error2_velocity_archive[i] = controller.error_velocity[1];

    controller.r[0] = controller.controller_u[0];  // bounded reference input of the first reference model
    controller.r[1] = controller.controller_u[2];  // bounded reference input of the second reference model
    controller.x1[0] = controller.controller_u[4]; // state variable of the first reference model
    controller.x1[1] = controller.controller_u[5];
    controller.x2[0] = controller.controller_u[6]; // state variable of the second reference model
    controller.x2[1] = controller.controller_u[7];

    for (int j = 0; j < n_joint; j++){
        controller.x[j] = controller.x1[j]; // state variable
        controller.x[j + 2] = controller.x2[j];
    }

    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint * 2; k++){
            controller.z[k][j] = controller.x[k]; // defined in Eq. 5
        }
        controller.z[4][j] = controller.r[j];
    }

    for (int j = 0; j < n_joint; j++){
        controller.error_state[0][j] = controller.error[j]; // tracking error of reference model
        controller.error_state[1][j] = controller.error_velocity[j];
    }

    // membership function of q1
    // first membership function of q1, V11, implies that q1 is small
    if (controller.x[0] <= 0 && controller.x[0] >= -PI / 2){
        membership.q1[0] = -controller.x[0] / (PI / 2);
    } else{
        membership.q1[0] = 0;
    }

    // second membership function of q1, V12, implies that q1 is medium
    if (controller.x[0] <= 0 && controller.x[0] >= -PI / 2){
        membership.q1[1] = (controller.x[0] + PI / 2) / (PI / 2);
    } else if (controller.x[0] <= PI && controller.x[0] >= 0){
        membership.q1[1] = (PI / 2 - controller.x[0]) / (PI / 2);
    } else{
        membership.q1[1] = 0;
    }

    // third membership function of q1, V13, implies that q1 is large
    if (controller.x[0] <= PI && controller.x[0] >= 0){
        membership.q1[2] = controller.x[0] / (PI / 2);
    } else{
        membership.q1[2] = 0;
    }

    // membership function of q2
    // first membership function of q2, V21, implies that q2 is small
    if (controller.x[2] <= 0 && controller.x[2] >= -PI / 2){
        membership.q2[0] = -controller.x[2] / (PI / 2);
    } else{
        membership.q2[0] = 0;
    }

    // second membership function of q2, V22, implies that q2 is medium
    if (controller.x[2] <= 0 && controller.x[2] >= -PI / 2){
        membership.q2[1] = (controller.x[2] + PI / 2) / (PI / 2);
    } else if (controller.x[2] <= PI && controller.x[2] >= 0){
        membership.q2[1] = (PI / 2 - controller.x[2]) / (PI / 2);
    } else{
        membership.q2[1] = 0;
    }

    // third membership function of q2, V23, implies that q2 is large
    if (controller.x[2] <= PI && controller.x[2] >= 0){
        membership.q2[2] = controller.x[2] / (PI / 2);
    } else{
        membership.q2[2] = 0;
    }

    // for (int j = 0; j < 3; j++) {
    //     printf("membership.q1[%d] = %lf\n", j, membership.q1[j]);
    // }

    double sum_q1 = 0, sum_q2 = 0;
    for (int j = 0; j < 3; j++){
        sum_q1 += membership.q1[j];
        sum_q2 += membership.q2[j];
    }

    for (int j = 0; j < 3; j++){
        controller.xi[j][0] = membership.q1[j] / sum_q1; // normalized firing strength vector of first reference model, defined in Eq. 6
        controller.xi[j][1] = membership.q2[j] / sum_q2; // normalized firing strength vector of second reference model, defined in Eq. 6
    }

    // calculate product of integral term's gamma, bc, solution of Lyapunov equation P, and error
    double gamma2_bc_P_error[n_joint];
    for (int j = 0; j < n_joint; j++){
        double sum = 0.0;
        for (int k = 0; k < n_joint; k++){
            sum += controller.P[1][k][0] * controller.error_state[k][j];
        }
        gamma2_bc_P_error[j] = controller.gamma[1][j] * sum;
        // printf("gamma2_bc_P_error[%d] = %lf\n", j, gamma2_bc_P_error[j]);
    }

    // derivative of integral term of PI update law, Eq. 9
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                controller.theta_derivative[k][g][j] = gamma2_bc_P_error[j] * controller.xi[k][j] * controller.z[g][j];
            }
        }
    }

    // for (int j = 0; j < n_joint; j++) {
    //     printf("controller.theta_derivative[%d]:\n", j);
    //     for (int k = 0; k < 3; k++) {
    //         for (int g = 0; g < n_joint * 2 + 1; g++) {
    //             printf("%8.4f ", controller.theta_derivative[k][g][j]);
    //             if ((g+1) % 5 == 0) {
    //                 printf("\n");
    //             }
    //         }
    //     }
    // }

    // handle non-numeric floating point value
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                if (isnan(controller.theta_derivative[k][g][j])) {
                    controller.theta_derivative[k][g][j] = 0;
                }
            }
        }
    }

    // integral term of parameter update PI law
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                controller.theta[k][g][j] += controller.theta_derivative[k][g][j] * Ts;
            }
        }
    }

    // calculate product of proportionality term's gamma, bc, solution of Lyapunov equation P, and error
    double gamma1_bc_P_error[n_joint];
    for (int j = 0; j < n_joint; j++){
        double sum = 0.0;
        for (int k = 0; k < n_joint; k++){
            sum += controller.P[1][k][1] * controller.error_state[k][j];
        }
        gamma1_bc_P_error[j] = controller.gamma[0][j] * sum;
        // printf("gamma1_bc_P_error[%d] = %lf\n", j, gamma1_bc_P_error[j]);
    }

    // proportionality term of PI update law, Eq. 9
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                controller.phi[k][g][j] = gamma1_bc_P_error[j] * controller.xi[k][j] * controller.z[g][j];
            }
        }
    }

    // fuzzy controller parameter update PI law, Eq. 37
    if (i == 0){
        for (int j = 0; j < n_joint; j++){
            for (int k = 0; k < 3; k++){
                for (int g = 0; g < n_joint * 2 +1; g++){
                    controller.THETA[k][g][j] = controller.phi[k][g][j];
                }
            }
        }
    } else{
        for (int j = 0; j < n_joint; j++){
            for (int k = 0; k < 3; k++){
                for (int g = 0; g < n_joint * 2 +1; g++){
                    controller.THETA[k][g][j] = controller.phi[k][g][j] + controller.theta[k][g][j];
                }
            }
        }
    }

    // handle of fuzzy controller parameter out-of-range
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                if (controller.THETA[k][g][j] >= 60){
                    controller.THETA[k][g][j] = 60;
                } else if (controller.THETA[k][g][j] <= -60){
                    controller.THETA[k][g][j] = -60;
                }
            }
        }
    }

    // handle non-numeric floating point value of fuzzy controller parameter
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < 3; k++){
            for (int g = 0; g < n_joint * 2 +1; g++){
                if (isnan(controller.THETA[k][g][j])) {
                    controller.THETA[k][g][j] = 0;
                }
            }
        }
    }

    // calculate product of normalized firing strength, fuzzy controller parameter
    double xi_THETA[n_joint * 2 +1][n_joint];
    for (int j = 0; j < n_joint; j++){
        for (int g = 0; g < n_joint * 2 +1; g++){
            double sum = 0.0;
            for (int k = 0; k < 3; k++){
                sum += controller.xi[k][j] * controller.THETA[k][g][j];
            }
            xi_THETA[g][j]= sum;
        }
    }

    // fuzzy logic system obtained through the combination of single-point fuzzification, product inference, and center-weighted defuzzification
    for (int j = 0; j < n_joint; j++){
        double sum = 0.0;
        for (int g = 0; g < n_joint * 2 +1; g++){
            sum += xi_THETA[g][j] * controller.z[g][j];
        }
        torque.torque_fuzzy[j] = sum; // output of fuzzy controller
        // printf("torque.torque_fuzzy[%d] = %lf\n", j, torque.torque_fuzzy[j]);
    }

    // handle non-numeric floating point value of output of fuzzy controller
    if (isnan(torque.torque_fuzzy[0])) {
        torque.torque_fuzzy[0] = controller.controller_u[8];
    }
    if (isnan(torque.torque_fuzzy[1])) {
        torque.torque_fuzzy[1] = controller.controller_u[9];
    }

    controller.Beta[0][0] = 0.8 * fabs(controller.x[3]); // known function
    controller.Beta[1][1] = 6.5 * fabs(controller.x[3]);

    // calculate product of error, solution of Lyapunov equation P, and bc
    double error_P_bc[n_joint];
    for (int j = 0; j < n_joint; j++){
        double sum = 0.0;
        for (int k = 0; k < n_joint; k++){
            sum += controller.error_state[k][j] * controller.P[k][1][j];
        }
        error_P_bc[j] = sum;
    }

    // calculate product of error, solution of Lyapunov equation P
    double error_P[n_joint][n_joint];
    for (int j = 0; j < n_joint; j++){
        for (int k = 0; k < n_joint; k++){
            double sum = 0.0;
            for (int g = 0; g < n_joint; g++){
                sum += controller.error_state[g][j] * controller.P[g][k][j];
            }
            error_P[k][j] = sum;
            // printf("error_P[%d][%d] = %lf\n", k, j, error_P[k][j]);
        }
    }

    // calculate product of error, solution of Lyapunov equation P and error
    double error_P_error[n_joint];
    for (int j = 0; j < n_joint; j++){
        double sum = 0.0;
        for (int k = 0; k < n_joint; k++){
            sum += error_P[k][j] * controller.error_state[k][j];
        }
        error_P_error[j] = sum;
    }

    // part of sliding mode control term
    double temp11, temp12, temp13, temp21, temp22, temp23;
    temp11 = (1 / controller.b[0][0][0] * (controller.b[0][1][1] * fabs(controller.controller_u[9])+ controller.d_upper[0])) * (1 + (controller.delta[0] / controller.beta[0]));
    temp12 = (error_P_bc[0] / (fabs(error_P_bc[0]) + controller.delta[0]));
    temp13 = controller.Beta[0][0] / (2 * pow(controller.b[0][0][0], 2)) * error_P_error[0];
    temp21 = (1 / controller.b[1][1][0] * (controller.b[1][0][1] * fabs(controller.controller_u[8]) + controller.d_upper[1])) * (1 + (controller.delta[1] / controller.beta[1]));
    temp22 = (error_P_bc[1] / (fabs(error_P_bc[1]) + controller.delta[1]));
    temp23 = controller.Beta[1][1] / (2 * pow(controller.b[1][1][0], 2)) * error_P_error[1];

    // additional control term used to overcome uncertainties
    torque.torque_si[0] =  temp11 * temp12 + temp13;
    torque.torque_si[1] =  temp21 * temp22 + temp23;

    // control input torque, Eq. 11
    for (int j = 0; j < n_joint; j++){
        torque.torque[j] = torque.torque_fuzzy[j] +  torque.torque_si[j];
        // printf("torque.torque[%d] = %lf\n", j, torque.torque[j]);
    }

    // handle of control input torque out-of-range
    if (torque.torque[0] >= 45){
        torque.torque[0] = controller.controller_u[8];
    } else if (torque.torque[0] <= -45){
        torque.torque[0] = controller.controller_u[8];
    }
    if (torque.torque[1] >= 45){
        torque.torque[1] = controller.controller_u[9];
    } else if (torque.torque[1] <= -45){
        torque.torque[1] = controller.controller_u[9];
    }

    // handle non-numeric floating point value of control input torque
    if (isnan(torque.torque[0])) {
        torque.torque[0] = controller.controller_u[8];
    }
    if (isnan(torque.torque[1])) {
        torque.torque[1] = controller.controller_u[9];
    }
    controller.controller_out[0] = torque.torque[0];
    controller.controller_out[1] = torque.torque[1];
    archive.torque1_archive[i] = torque.torque[0];
    archive.torque2_archive[i] = torque.torque[1];
}

struct _plant{
    double plant_u[2];
    double plant_out[4];
} plant;

void PLANT_init(){
    system_state.q[0] = 0.1;
    system_state.dq[0] = 2.0;
    system_state.q[1] = 0.1;
    system_state.dq[1] = 2.5;
    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.q[1];
    plant.plant_out[3] = system_state.dq[1];
    l = 1.0;     // length of linl
    m[0] = 1.0;  // mass of link 1
    m[1] = 2.0;  // mass of link 2
    kv[0] = 0.3; // coefficient of viscous friction
    kv[1] = 0.5;
    kc[0] = 0.2; // coefficient of coulomb friction
    kc[1] = 0.5;
}

double PLANT_realize(int i){
    plant.plant_u[0] = controller.controller_out[0];
    plant.plant_u[1] = controller.controller_out[1];

    // inertia matrix for manipulator dynamics equation
    dynamics.D[0][0] = 1.0/3.0 * m[0] * pow(l, 2) + 4.0/3.0 * m[1] * pow(l, 2) + m[1] * pow(l, 2) * cos(controller.x[2]);
    dynamics.D[0][1] = 1.0/3.0 * m[1] * pow(l, 2) + 1.0/2.0 * m[1] * pow(l, 2) * cos(controller.x[2]);
    dynamics.D[1][0] = dynamics.D[0][1];
    dynamics.D[1][1] = 1.0/3.0 * m[1] * pow(l, 2);

    // Coriolis and centrifugal matrix for manipulator dynamics equation 
    dynamics.C[0][0] = -m[1] * pow(l, 2) * sin(controller.x[2]) * controller.x[3];
    dynamics.C[0][1] = -1.0/2.0 * m[1] * pow(l, 2) * sin(controller.x[2]) * controller.x[3];
    dynamics.C[1][0] = 1.0/2.0 * m[1] * pow(l, 2) * sin(controller.x[2]) * controller.x[1];
    dynamics.C[1][1] = 0.0;

    // gravitational torque matrix for manipulator dynamics equation
    dynamics.G[0] = 1.0/2.0 * m[0] * g * l * cos(controller.x[0]) + 1.0/2.0 * m[1] * g * l * cos(controller.x[0] + controller.x[2]) + m[1] * g * l * cos(controller.x[0]);
    dynamics.G[1] = 1.0/2.0 * m[1] * g * l * cos(controller.x[0] + controller.x[2]);

    // viscous friction torque
    dynamics.Fv[0] = kv[0] * controller.x[1];
    dynamics.Fv[1] = kv[1] * controller.x[3];

    // coulomb friction torque
    dynamics.Fc[0] = kc[0] * sign(controller.x[1]);
    dynamics.Fc[1] = kc[1] * sign(controller.x[3]);

    // for (int j = 0; j < 2; j++) {
    //     printf("dynamics.Fc[%d]: %lf\n", j, dynamics.Fc[j]);
    // }

    double inv_D[n_joint][n_joint], C_dq[n_joint], torque_Cdq_G_Fv_Fc[n_joint];
    inv_matrix(inv_D, dynamics.D, n_joint); // calculate inverse of inertia matrix for manipulator dynamics equation
    matrix_Multi((double *)C_dq, (double *)dynamics.C, (double *)system_state.dq, 2, 2, 1); // calculate C multiplied by actual angular velocity

    for (int j = 0; j < n_joint; j++){
        torque_Cdq_G_Fv_Fc[j] = torque.torque[j] - C_dq[j] - dynamics.G[j] - dynamics.Fv[j] - dynamics.Fc[j];
    }

    // actual manipulator dynamics system
    matrix_Multi((double *)system_state.ddq, (double *)inv_D, (double *)torque_Cdq_G_Fv_Fc, 2, 2, 1);

    system_state.dq[0] = system_state.dq[0] + system_state.ddq[0] * Ts;
    system_state.dq[1] = system_state.dq[1] + system_state.ddq[1] * Ts;
    system_state.q[0] = system_state.q[0] + system_state.dq[0] * Ts;
    system_state.q[1] = system_state.q[1] + system_state.dq[1] * Ts;

    plant.plant_out[0] = system_state.q[0];
    plant.plant_out[1] = system_state.dq[0];
    plant.plant_out[2] = system_state.q[1];
    plant.plant_out[3] = system_state.dq[1];
}

void saveArchiveToTxt(double *archive, int size, const char *filename){

    FILE *file = fopen(filename, "w+");

    if (file == NULL){
        perror("Failed to open file");
        exit(1);
    }
    else{
        for (int i = 0; i < size; i++){
            fprintf(file, "%lf\n", archive[i]);
        }
        fclose(file);
        printf("Saved to file %s\n", filename);
    }
}

void saveArchive(){

    saveArchiveToTxt(q1_desired.y, ARRAY_SIZE, "../report/qd1.txt");
    saveArchiveToTxt(dq1_desired.y, ARRAY_SIZE, "../report/dqd1.txt");
    saveArchiveToTxt(archive.q1_archive, ARRAY_SIZE, "../report/q1.txt");
    saveArchiveToTxt(archive.dq1_archive, ARRAY_SIZE, "../report/dq1.txt");
    saveArchiveToTxt(q2_desired.y, ARRAY_SIZE, "../report/qd2.txt");
    saveArchiveToTxt(dq2_desired.y, ARRAY_SIZE, "../report/dqd2.txt");
    saveArchiveToTxt(archive.q2_archive, ARRAY_SIZE, "../report/q2.txt");
    saveArchiveToTxt(archive.dq2_archive, ARRAY_SIZE, "../report/dq2.txt");
    saveArchiveToTxt(archive.error1_archive, ARRAY_SIZE, "../report/error1.txt");
    saveArchiveToTxt(archive.error1_velocity_archive, ARRAY_SIZE, "../report/error1_velocity.txt");
    saveArchiveToTxt(archive.error2_archive, ARRAY_SIZE, "../report/error2.txt");
    saveArchiveToTxt(archive.error2_velocity_archive, ARRAY_SIZE, "../report/error2_velocity.txt");
    saveArchiveToTxt(archive.torque1_archive, ARRAY_SIZE, "../report/torque1.txt");
    saveArchiveToTxt(archive.torque2_archive, ARRAY_SIZE, "../report/torque2.txt");
}

int main(){

    SystemInput(&q1_desired, &qd1_1, &qd1_2, &dq1_desired, &dqd1_1, &dqd1_2, &q2_desired, &qd2_1, &qd2_2, &dq2_desired, &dqd2_1, &dqd2_2, Ts, t0, t1);
    CONTROLLER_init(); // initialize controller parameter
    PLANT_init();      // initialize plant parameter

    for (int i = 0; i < ARRAY_SIZE; i++){
    // for (int i = 0; i < 5; i++){
        double time = i * Ts + t0;
        printf("time at step %d: %f\n", i, time);
        CONTROLLER_realize(i);
        PLANT_realize(i);
    }

    saveArchive();

    return 0;
}
