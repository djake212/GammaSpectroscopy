import sympy as sp
import numpy as np

def error_propagation(f, var_values, var_errors):
    """
    f: function with sympy symbols as input variables
    var_values: dictionary with keys as the sympy symbols and values as the variable values
    var_errors: dictionary with keys as the sympy symbols and values as the variable errors
    """

    # Calculate the partial derivatives of the function for each variable
    partial_derivatives = [sp.diff(f, var) for var in var_values.keys()]

    # Calculate the propagated error using partial derivatives
    sigma_f_squared = sum((partial_derivative * var_errors[var])**2 for partial_derivative, var in zip(partial_derivatives, var_values.keys()))

    # Substitute the variable values
    sigma_f_val = sp.sqrt(sigma_f_squared.subs(var_values))

    return sigma_f_val

def usrprompt():
	usrinput = input("Enter a date in the form DDMmmYYYY: ")
	return usrinput
	
#function to parse input string and return list containing year, month number, and day of the month
def dateparse(datestr):
	strlst = []
	for char in datestr:
		strlst.append(char)
	datelst = []
	daynum = int(strlst[0] + strlst[1])
	datelst.append(daynum)
	monthdict = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
	monthnum = monthdict[strlst[2]+strlst[3]+strlst[4]]
	datelst.append(monthnum)
	yearnum = int(strlst[5]+strlst[6]+strlst[7]+strlst[8])
	datelst.append(yearnum)

	return datelst

#function to compute julian day number using returned list
def juliancalc(datelst):
	D = datelst[0]
	M = datelst[1]
	Y = datelst[2]
	if M <= 2:
		Y -= 1
		M += 12
	A = int(Y/100)
	B = 2 - A + int(A/4)
	JD = int(365.25*(Y+4716)) + int(30.6001 * (M + 1)) + D + B - 1524.5
	return JD

def values_agree(value1, error1, value2, error2):
    min_value1 = value1 - error1
    max_value1 = value1 + error1
    min_value2 = value2 - error2
    max_value2 = value2 + error2

    # Check if the error ranges overlap
    if max_value1 >= min_value2 and max_value2 >= min_value1:
        return True
    else:
        return False

# Define the var and functions
A0, Thalf, time = sp.symbols('A0 Thalf time')
Rk, Tk, Tn, Nk, Nn = sp.symbols('Rk Tk Tn Nk Nn')
findcurrentactivity_sympy = A0 * (3.7 * 10**4) * sp.exp((-time * sp.log(2)) / Thalf) / (3.7 * 10**4)
activity_sympy = Rk * (Tk / Tn) * (Nn / Nk)

#create dictionary for the samples
#enter days elapsed since april 12th if new data is taken
elapsed = 5

samples = {
    "Na22": {
        "halflife": 2.6*365,
        "halflife_error": 0.1*365,
        "activity": 1,
        "activity_error": 0.2,
        "date": 160.0 + elapsed,
        "date_error": 1,
    },
    "Ba133": {
        "halflife": 10.6*365,
        "halflife_error": 0.1*365,
        "activity": 1,
        "activity_error": 0.2,
        "date": 84.0 + elapsed,
        "date_error": 1,
    },
    "Co60": {
        "halflife": 5.27*365,
        "halflife_error": 0.01*365,
        "activity": 1,
        "activity_error": 0.2,
        "date": 1490.0 + elapsed,
        "date_error": 15,
    },
    "Mn54": {
        "halflife": 312,
        "halflife_error": 1,
        "activity": 1,
        "activity_error": 0.2,
        "date": 1490 + elapsed,
        "date_error": 15,
    },
    "Co57": {
        "halflife": 271.74,
        "halflife_error": 0.01,
        "activity": 1,
        "activity_error": 0.2,
        "date": 160 + elapsed,
        "date_error": 1,
    },
    "Cs137": {
        "halflife": 30.08*365,
        "halflife_error": 0.01*365,
        "activity": 0.25,
        "activity_error": 0.05,
        "date": 162 + elapsed,
        "date_error": 1,
    },
    "Cd109": {
        "halflife": 461.4,
        "halflife_error": 0.1,
        "activity": 1,
        "activity_error": 0.2,
        "date": 117 + elapsed,
        "date_error": 1,
    },
    "Eu152": {
        "halflife": 13.52*365,
        "halflife_error": 0.01*365,
        "activity": 0.5,
        "activity_error": 0.1,
        "date": 217 + elapsed,
        "date_error": 1,
    },
    "Zn65": {
        "halflife": 243.93,
        "halflife_error": 0.01,
        "activity": 1,
        "activity_error": 0.2,
        "date": 217.0 + elapsed,
        "date_error": 1
    },
    "Unknown": {
        "halflife": 0,
        "halflife_error": 0,
        "activity": 0.5,
        "activity_error": 0.1,
        "activity2": 1,
        "activity_error2": 0.2,
        "date": 1490.0 - 5 + elapsed,
        "date_error": 15
    }
}

#known integral
knownintegral= input("Paste output for known sample gaussian: ")
start_index = knownintegral.index("integral value")
integral_substring = knownintegral[start_index:]
integral_parts = integral_substring.split("(error:")
integral_value_str = integral_parts[0].split("=")[1].strip()
integral_error_str = integral_parts[1].split(")")[0].strip("+-")
integral_value = float(integral_value_str)
integral_error = float(integral_error_str)

#unknownintegral
unknownintegral= input("Paste output for unknown sample gaussian: ")
start_index2 = unknownintegral.index("integral value")
integral_substring2 = unknownintegral[start_index2:]
integral_parts2 = integral_substring2.split("(error:")
integral_value_str2 = integral_parts2[0].split("=")[1].strip()
integral_error_str2 = integral_parts2[1].split(")")[0].strip("+-")
integral_value2 = float(integral_value_str2)
integral_error2 = float(integral_error_str2)

#enter known sample to compare
knownsample = input('Enter which known sample to use: ').strip().capitalize()
timeknown = float(input('Enter the Live time for the known sample: ').strip())
timeunknown = float(input('Enter the Live time for the unknown sample: ').strip())

var_values = {A0: samples[knownsample]['activity'], Thalf: samples[knownsample]['halflife'], time: samples[knownsample]['date']}
var_errors = {A0: samples[knownsample]['activity_error'], Thalf: samples[knownsample]['halflife_error'], time: samples[knownsample]['date_error']}
Rkval = findcurrentactivity_sympy.subs([(A0,var_values.get(A0)),(Thalf,var_values.get(Thalf)),(time,var_values.get(time))])
print()
print(f'Reactivity of known sample activity currently in microCi: {Rkval}')

# Calculate propagated error for findcurrentactivity function
error_findcurrentactivity = error_propagation(findcurrentactivity_sympy, var_values, var_errors)
error_findcurrentactivity = error_findcurrentactivity.evalf()
print(f"Error in  known sample activitycurrently: {error_findcurrentactivity} and as percent: {100*error_findcurrentactivity/Rkval:.2f}%")
print()

#Enter variables for current sample

var_values2 = {Rk: Rkval, Tk: timeknown, Tn: timeunknown, Nk: integral_value, Nn: integral_value2}
var_errors2 = {Rk: error_findcurrentactivity, Tk: 1, Tn: 1, Nk: integral_error, Nn: integral_error2}

# Calculate propagated error for activity function
error_activity = error_propagation(activity_sympy, var_values2, var_errors2)
error_activity = error_activity.evalf()
Rkval2 = activity_sympy.subs([(Rk, var_values2.get(Rk)), (Tk, var_values2.get(Tk)), (Tn, var_values2.get(Tn)), (Nk, var_values2.get(Nk)), (Nn, var_values2.get(Nn))])
print(f'Experimental activity of unknown in microCi: {Rkval2}')
print(f"Error in activity of unknown: {error_activity} and as percent: {100*error_activity/Rkval2:.2f}%")
print()

#compare to theory for unknown change the A0 to either activity or activity2 to test both activities listed
var_values3 = {A0: samples['Unknown']['activity'], Thalf: samples[knownsample]['halflife'], time: samples['Unknown']['date']}
var_errors3 = {A0: samples['Unknown']['activity_error'], Thalf: samples[knownsample]['halflife_error'], time: samples['Unknown']['date_error']}
Rkval3 = findcurrentactivity_sympy.subs([(A0,var_values3.get(A0)),(Thalf,var_values3.get(Thalf)),(time,var_values3.get(time))])
print(f'Theory Reactivity of current in microCi: {Rkval3}')

error_findcurrentactivity2 = error_propagation(findcurrentactivity_sympy, var_values3, var_errors3)
error_findcurrentactivity2 = error_findcurrentactivity2.evalf()

print(f"Error in theory of inital activity: {error_findcurrentactivity2} and as percent: {100*error_findcurrentactivity2/Rkval3:.2f}%")
print(f'the % difference between the experimental and expected value is {100*(Rkval3-Rkval2)/Rkval3:.2f}%')
valueagree = values_agree(Rkval3,error_findcurrentactivity2,Rkval2,error_activity)
if valueagree:
	print('The values are within experimental error!')
else:
	print('These values are not within experimental error, something went wrong somewhere :(')
print()