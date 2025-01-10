import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
import matpopmod as mpm

# Exponential model parameters:
exponential_rate = 1.4538071538944615  # exponential rate of growth: value of k in p(t) = p_0 * e^(kt)

# Logistic model parameters:
initial_population = 100
birth_rate_slope = 0.0001  # rate of change of birth rate with respect to population (negative)
birth_rate_max = 20  # birth rate at population = 0
death_rate = 19  # constant death rate of population

# Matrix model parameters:
old_transition_matrix = np.array([[0, 1, 1.9, 1.7, 0],
                                  [0.61, 0, 0, 0, 0],
                                  [0, 0.33, 0, 0, 0],
                                  [0, 0, 0.21, 0, 0],
                                  [0, 0, 0, 0.13, 0.08]])  # transition matrix
new_fish_tech = np.array([[0, 0, 0.1, 0, 0],
                          [0, 0, 0, 0, 0],
                          [0, 0.04, 0, 0, 0],
                          [0, 0, 0.02, 0, 0],
                          [0, 0, 0, 0, 0]])  # changes to transition_matrix due to new fishing technology
new_trans_matrix = np.subtract(old_transition_matrix, new_fish_tech)  # new_fish_tech applied to old_trans_matrix


def change_survival(tram, survivals):
    f_changes = [0, 0, 0, 0, 0]
    for i in range(len(tram[0]) - 1):
        f_changes[i] = tram[0][i] - tram[0][i] * (1 + survivals[i] / tram[i+1][i])
    f_changes[len(tram[0]) - 1] = tram[0][len(tram[0]) - 1] * (1 + survivals[len(tram[0]) - 1] / tram[len(tram[0]) - 1][len(tram[0]) - 1])
    matrix = [f_changes,
        [survivals[0], 0, 0, 0, 0],
        [0, -survivals[1], 0, 0, 0],
        [0, 0, -survivals[2], 0, 0],
        [0, 0, 0, -survivals[3], -survivals[4]]]  # stable fishing laws
    return np.array(matrix)


stable_fishing_changes = change_survival(old_transition_matrix, [0, -0.011, -0.011, -0.0115547474741, 0])
stable_trans_matrix = np.subtract(old_transition_matrix, stable_fishing_changes)

print(stable_trans_matrix)

pop_vec = [289, 211, 120, 76, 51]  # initial population vector
# pop_vec = [407, 248, 79, 16, 2]


# Constants
initial_population_rate = (birth_rate_slope * initial_population) * ((birth_rate_max - death_rate) / birth_rate_slope -
                                                                     initial_population)  # determines initial rate of change of population for logistic model
p_0e = initial_population, exponential_rate  # initial value of "p" for exponential model
p_0l = initial_population, initial_population_rate  # initial value of "p" for logistic model

# General parameters:
time_interval = 5000  # number of time steps the model will run for
time_steps = 10000  # clarity of graph: higher number = more points plotted


def exponential(p, k):  # exponential model function
    result = [p]
    for i in range(time_steps - 1):
        t = time_interval * (i / time_steps)
        pop_i = p * (np.e ** (np.log(k) * t))
        result.append(pop_i)
    return result


def logistic(p, t):  # logistic model function
    population, population_rate = p
    population_rate = (birth_rate_slope * population) * ((birth_rate_max - death_rate) / birth_rate_slope - population)
    population_rate_rate = population_rate * (birth_rate_max - death_rate) - (
            2 * population * population_rate * birth_rate_slope)
    return population_rate, population_rate_rate


# Exponential Model
# p_e = exponential(initial_population, exponential_rate)  # run exponential function
# plt.plot(time, p_e)  # plot exponential model data

# Logistic Model
# time = np.linspace(0, time_interval, time_steps)  # format time for odeint()
# p_l = odeint(logistic, p_0l, time)  # run logistic model
# plt.plot(time, p_l[:, 0])  # plot logistic model data

# Matrix model
# old transition matrix
# old_tram = mpm.MPM(A=old_transition_matrix)
# traj = old_tram.trajectory(pop_vec, t_max=time_interval)
# plot1 = mpm.plot.trajectory(traj, second_order=False)
# plot1.set(xlabel="Time-steps", ylabel="Cod Population")
# old_trajs = old_tram.stochastic_trajectories(pop_vec, t_max=time_interval, reps=1000)
# mpm.plot.multiple_trajectories(old_trajs)

# new transition matrix
# new_tram = mpm.MPM(A=new_trans_matrix)
# traj = new_tram.trajectory(pop_vec, t_max=time_interval)
# plot1 = mpm.plot.trajectory(traj, second_order=False)
# plot1.set(xlabel="Time-steps", ylabel="Cod Population")
# new_trajs = new_tram.stochastic_trajectories(pop_vec, t_max=time_interval, reps=1000)
# mpm.plot.multiple_trajectories(new_trajs)

# stable transition matrix
stable_tram = mpm.MPM(A=stable_trans_matrix)
traj = stable_tram.trajectory(pop_vec, t_max=time_interval)
plot1 = mpm.plot.trajectory(traj, second_order=False)
plot1.set(xlabel="Time-steps", ylabel="Cod Population")
stable_trajs = stable_tram.stochastic_trajectories(pop_vec, t_max=time_interval, reps=1000)
stab_plot = mpm.plot.multiple_trajectories(stable_trajs)
stab_plot.set(xlabel="Time-steps", ylabel="Cod Population")

# print("\nMatrix model:\n"
#       "Stable stage distribution (A):", stable_tram.w,
#       "\nGrowth rate (main):", stable_tram.lmbd)

# Show plots
plt.show()
mpm.plot.show()

# print("Logistic model:\n"
#       "k =", birth_rate_slope,
#       "\nM =", (birth_rate_max - death_rate) / birth_rate_slope)
