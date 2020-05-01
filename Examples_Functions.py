"""
Michael Hollas - 2020-05-01
------------------------------------------------------------------------------------------------------------------------
How To use functions
------------------------------------------------------------------------------------------------------------------------

functions are created by defining them (hence def)

def functionname(input variables):
"""
def simpleAddition(a,b):
    #do the thing!
    return a + b # will return this value to be the output of the function - don't need to have a return





print(simpleAddition(1,2))

# or

addValue = simpleAddition(1,2)

print(addValue)


"""
input/arguments with keywords! 
"""

print(simpleAddition(a=2,b=5))

#order does not matter

print(simpleAddition(b=2,a=5))


"""
default values!
"""

def my_function(firstname = "john",lastname="smith"):
    print("{0} {1}".format(firstname,lastname))

my_function()



"""
variables inside a function usually cannot be accessed outside the function 
"""
def simplemultiplaction(a,b):
    x = 12
    return a*b


#print(x) #<--- can't do!!!'


# functions can be called from another file!!!! :O

"""
first need to import the function
"""

from Examples_CSV_file_reader import csv_to_dict as cd

# or

from Examples_CSV_file_reader import * # this will import all functions from this file


from Examples_CSV_file_reader import csv_to_dict as cd # can rename the function with the 'as' keyword



data = csv_to_dict("C:\Data\\test.csv")

data2 = cd("C:\Data\\test.csv")# this is the same function as above!



