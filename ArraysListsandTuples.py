"""
Michael Hollas - 2020-05-01
------------------------------------------------------------------------------------------------------------------------
Arrays, Lists and Tuples
------------------------------------------------------------------------------------------------------------------------
"""




"""
Arrays don't Exist! don't blame me! you can use lists or use the Numpy library
"""


"""
Lists
"""

numbers = [1,2,3,4,5,6]  # lists are defined using square brackets

list_of_strings = ["hello","world","my","name","is", "Mike"]

mixed_list = ["HCD", 2344.5324,"MCRRTP"] # you can have mixed data types in lists!!!!! suck on that c#!!

#you can access an element via indexing (lists and tuples are 0 indexed! - the first one is 0)

print(mixed_list[0]) # this will print the first element of the list

#fancy indexing!!!

print(numbers[2:]) # yields a new list including the 3rd value to the end of the list

print(numbers[:3]) # yields first 4 elements in the list

print(numbers[2:3]) #

print(numbers[-1]) # gives the last element of the list!!!


#changing values in a list

mixed_list[0] = "HeLOO"
print(mixed_list)



"""
MOST IMPORTANT THING ABOUT LISTS!

They do not have a defined size when created - you can APPEND to them
"""

mixed_list.append("new thing")
print(mixed_list)



"""
looping through a list
"""

for i in range(len(mixed_list)): # indexed based version
    print(mixed_list[i])


#or

for thing in mixed_list: # for each based version
    print(thing)

"""
these do the same 
"""


