import numpy
import pandas
import matplotlib.pyplot as plt
from tqdm.contrib.itertools import product
from sklearn.metrics import mean_squared_error


def sird_forecast(beta, gamma, mu, step, nb_jours, init_conditions):
    """
    Simulate the SIRD model.
    
    Parameters:
        beta, gamma, mu: SIRD model parameters.
        step: Time step for the simulation.
        nb_jours: Total duration (in days).
        init_conditions: Dictionary with keys "S", "I", "R", "D".
        
    Returns:
        time: List of time points.
        S, I, R, D: Arrays of population fractions.
    """
    time = [0]
    S = [init_conditions["S"]]
    I = [init_conditions["I"]]
    R = [init_conditions["R"]]
    D = [init_conditions["D"]]
    
    nb_points = int(nb_jours / step)
    
    for _ in range(1, nb_points):
        # current values
        s_current = S[-1]
        i_current = I[-1]
        r_current = R[-1]
        d_current = D[-1]
        
        # Compute derivatives
        ds = -beta * s_current * i_current
        di = beta * s_current * i_current - gamma * i_current - mu * i_current
        dr = gamma * i_current
        dd = mu * i_current
        
        # Euler integration step
        s_next = s_current + ds * step
        i_next = i_current + di * step
        r_next = r_current + dr * step
        d_next = d_current + dd * step
        
        time.append(time[-1] + step)
        S.append(s_next)
        I.append(i_next)
        R.append(r_next)
        D.append(d_next)
    
    # Downsample for plotting (if needed)
    time = time[::1000]
    S = numpy.array(S[::1000])
    I = numpy.array(I[::1000])
    R = numpy.array(R[::1000])
    D = numpy.array(D[::1000])
    return time, S, I, R, D

def rmse(predictions, targets):
    return numpy.sqrt(mean_squared_error(targets, predictions))


def grid_search_sird(step, nb_jours, ground_truth, init_conditions):
    # Define parameter ranges (you may need to adjust these ranges)
    betas = numpy.linspace(0.1, 1, 20)
    gammas = numpy.linspace(0.1, 1, 4)
    mus = numpy.linspace(0.01, 0.5, 4)  # mortality is often smaller than the others


    #Best parameters found:
    #beta: 0.33684210526315794 gamma: 0.1 mu: 0.01
    best_beta, best_gamma, best_mu = None, None, None 
    best_rmse = float("inf")
    
    for beta, gamma, mu in product(betas, gammas, mus):
        time, S, I, R, D = sird_forecast(beta, gamma, mu, step, nb_jours, init_conditions)
        # Compute RMSE for each compartment. Ensure that ground_truth is aligned in time.
        rmse_S = rmse(S, ground_truth["Susceptibles"].values)
        rmse_I = rmse(I, ground_truth["Infectés"].values)
        rmse_R = rmse(R, ground_truth["Rétablis"].values)
        rmse_D = rmse(D, ground_truth["Décès"].values)
        actual_rmse = rmse_S + rmse_I + rmse_R + rmse_D
        
        if actual_rmse < best_rmse:
            best_rmse = actual_rmse
            best_beta, best_gamma, best_mu = beta, gamma, mu

    print("Best parameters found:")
    print("beta:", best_beta, "gamma:", best_gamma, "mu:", best_mu)
    
    time, best_S, best_I, best_R, best_D = sird_forecast(best_beta, best_gamma, best_mu, step, nb_jours, init_conditions)
    plot_sird(time, best_S, best_I, best_R, best_D, ground_truth)
    
    return time, best_S, best_I, best_R, best_D


def plot_sird(time, S, I, R, D, ground_truth):
    plt.figure(figsize=(15, 6))
    
    # Model predictions
    plt.plot(time, S, "b:", label='Susceptibles (predicted)')
    plt.plot(time, I, "r:", label='Infectés (predicted)')
    plt.plot(time, R, "g:", label='Rétablis (predicted)')
    plt.plot(time, D, "k:", label='Décès (predicted)')
    
    # Ground truth data (assuming the CSV has a column 'Jour' matching the sampled time points)
    plt.plot(ground_truth["Jour"], ground_truth["Susceptibles"], "b--", label='Susceptibles (real)')
    plt.plot(ground_truth["Jour"], ground_truth["Infectés"], "r--", label='Infectés (real)')
    plt.plot(ground_truth["Jour"], ground_truth["Rétablis"], "g--", label='Rétablis (real)')
    plt.plot(ground_truth["Jour"], ground_truth["Décès"], "k--", label='Décès (real)')
    
    plt.xlabel('Temps (Jours)')
    plt.ylabel('Population Fraction')
    plt.title('Modélisation SIRD et Données Réelles')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    step = 0.001    # time step in days
    nb_jours = 90   # simulation for 90 days

    # Load the ground truth data (ensure the CSV has 90 rows corresponding to days 0 to 89)
    ground_truth = pandas.read_csv("dev back\Lotka-Volterra\LotkaVolterra\sird_dataset.csv")
    
    # Define initial conditions.
    # For the epidemic model, you might want to adjust initial values if needed.
    init_conditions = {
        "S": ground_truth.loc[ground_truth["Jour"] == 0, "Susceptibles"].values[0],
        "I": ground_truth.loc[ground_truth["Jour"] == 0, "Infectés"].values[0],
        "R": 0,  # or the actual value if known
        "D": 0   # or the actual value if known
    }
    
    grid_search_sird(step, nb_jours, ground_truth, init_conditions)
