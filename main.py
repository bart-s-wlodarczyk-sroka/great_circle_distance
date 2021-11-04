import csv
import math
import pandas as pd
import numpy as np
import scipy.constants as con
from scipy.optimize import minimize

possible_targets = "/Users/Bart/Downloads/m13_max_pm_test.csv"
c_ra = 250.423475
c_dec = 36.46131944

def gcd(fun):
    delta_theta = np.arccos(np.sin(fun[1] + fun[3] * fun[4] / con.c) * np.sin(c_dec) + np.cos(fun[1] + fun[3] * fun[4] / con.c) * np.cos(c_dec) * np.cos(fun[0] + fun[2] * fun[4] / con.c - c_ra))
    return delta_theta

def neg_gcd(fun):
    delta_theta = (np.arccos(np.sin(fun[1] + fun[3] * fun[4] / con.c) * np.sin(c_dec) + np.cos(fun[1] + fun[3] * fun[4] / con.c) * np.cos(c_dec) * np.cos(fun[0] + fun[2] * fun[4] / con.c - c_ra))) * -1
    return delta_theta

source_ids = []
min_deltas = []
max_deltas = []

with open(possible_targets, "r") as csvfile:
    datareader = csv.reader(csvfile)
    next(datareader)
    count = 1
    for row in datareader:
        source_id = row[0]
        # COORDINATES IN DEG, ERRORS IN MAS
        ra_min = math.radians(float(row[1]) - (3.6e6 * float(row[2])))
        ra_max = math.radians(float(row[1]) + (3.6e6 * float(row[2])))
        dec_min = math.radians(float(row[3]) - (3.6e6 * float(row[4])))
        dec_max = math.radians(float(row[3]) + (3.6e6 * float(row[4])))
        # PM IN MAS/YR, PM ERROR IN MAS/YR
        pmra_min = math.radians(3.6e6 * (float(row[7]) - float(row[8]))) / 3.154e7
        pmra_max = math.radians(3.6e6 * (float(row[7]) + float(row[8]))) / 3.154e7
        pmdec_min = math.radians(3.6e6 * (float(row[9]) - float(row[10]))) / 3.154e7
        pmdec_max = math.radians(3.6e6 * (float(row[9]) + float(row[10]))) / 3.154e7
        # DISTANCE IN PC
        r_lo = float(row[21]) * 3.086e+16
        r_hi = float(row[22])* 3.086e+16

        ra_guess = math.radians(float(row[1]))
        dec_guess = math.radians(float(row[3]))
        pmra_guess = math.radians(3.6e6 * (float(row[7]))) / 3.154e7
        pmdec_guess = math.radians(3.6e6 * (float(row[9]))) / 3.154e7
        r_guess = float(row[20]) * 3.086e+16

        # print("ra_guess: {}".format(ra_guess))
        # print("dec_guess: {}".format(dec_guess))
        # print("pmra_guess: {}".format(pmra_guess))
        # print("pmdec_guess: {}".format(pmdec_guess))
        # print("r_guess: {}".format(r_guess))

        start_guess = [ra_guess, dec_guess, pmra_guess, pmdec_guess, r_guess]

        my_bounds = ((ra_min, ra_max), (dec_min, dec_max), (pmra_min, pmra_max), (pmdec_min, pmdec_max), (r_lo, r_hi))

        n_min = minimize(gcd, x0=start_guess, method='SLSQP', bounds=my_bounds)
        n_max = minimize(neg_gcd, x0=start_guess, method='SLSQP', bounds=my_bounds)

        print("Minimized RA: {}".format(math.degrees(n_min['x'][0])))
        print("Minimized DEC: {}".format(math.degrees(n_min['x'][1])))
        print("Minimum radial offset: {} rad".format(n_min["fun"]))
        print("Maximized RA: {}".format(math.degrees(n_max['x'][0])))
        print("Maximized DEC: {}".format(math.degrees(n_max['x'][1])))
        print("Maximum radial offset: {} rad".format(n_max["fun"] * -1))

        source_ids.append(source_id)
        min_deltas.append(n_max["fun"])
        max_deltas.append(n_min["fun"])

        print("Target {} completed.".format(count))
        count += 1

with open("/Users/Bart/Downloads/meti_error_output.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    headers = ["source_id", "min_delta_theta", "max_delta_theta"]
    writer.writerow(headers)
    data = list(zip(source_ids, min_deltas, max_deltas))
    for row in data:
        row = list(row)
        writer.writerow(row)
    print("All rows written to 'meti_error_output.csv'.")

# FOR EACH ROW IN THE CSV FILE,
# CREATE DATA STRUCTURE WITH SOURCE_ID
# APPEND MIN(DELTA_THETA) FOR EACH SOURCE_ID
# APPEND MAX(DELTA_THETA) FOR EACH SOURCE_ID
# CALCULATE N_MAX
# CALCULATE N_MIN
