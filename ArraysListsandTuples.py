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


#fancy indexing!!!




#changing values in a list

mixed_list[0] = "HeLOO"




"""
MOST IMPORTANT THING ABOUT LISTS!

They do not have a defined size when created - you can APPEND to them
"""

mixed_list.append("new thing")

mixed_list.insert(2,"i'm first now!!!!")






"""
looping through a list
"""

data_list = [["mz","int"],["rtwetw","ghdfhfh"],["djhkh","isdgg"],["jhhn","wercdf"]]


print(data_list[0][0])



with open("test.csv","w") as file:
    for i in range(len(data_list)):
        for j in range(len(data_list[i])):
            file.write(data_list[i][j] + ",")
        file.write("\n")






"""
for i in range(len(mixed_list)): # indexed based version
    print(mixed_list[i])


#or

for thing in mixed_list: # for each based version
    print(thing)
"""
"""
these do the same 
"""



"""
Tuples - mostly the same, but cannot be changed after creating 
"""



