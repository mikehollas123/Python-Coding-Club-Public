"""
Michael Hollas - 2020-04-24
------------------------------------------------------------------------------------------------------------------------

Here's and example of a csv to dictionary file reader.

the reader is contained in a function (I'll have another tutorial about these)

This will take the desired csv file, read it into memory and create a dictionary containing the data as columns.
If the first line of the file contains the column names, the user can signify this and the dictionaries keys
will be the names of these columns, otherwise the keys will be 'col 0' etc.

this is not the only way to do this - sometimes you may want data to be grouped by row - which will require
major modification to the how the code works - so should be created as a separate function

EVERYTHING IS STORED AS STRINGS - must be converted if you want to do any data manipulation
    How can we fix this? surely we can make something that will automatically determine the data type?!
        hint: look at try/except statements
------------------------------------------------------------------------------------------------------------------------
"""

"""
creates a function that can be called elsewhere. The input variables are contained in the parentheses. 

    Default values can be defined by defining the value of the variable (i.e. first_row_is_header=True). If no variable
    is given by the user, these values will be used.
"""


def csv_to_dict(file_location, first_row_is_header=True):  # creates the funtion

    with open(file_location,
              "r") as file:  # opens the file and creates a file object (don't worry if you don't understand objects yet)
        data = file.read()  # reads the whole file into memory as a single long string

        # lets split the string by each line in the file
        lines = data.split(
            "\n")  # sometimes need \n\r  ! varies from csv to csv. How can we account for this? (hint: regex)

        csv_dict = {}  # create an empty dictionary

        headings = lines[0].split(
            ",")  # split the first line to yield the headings (or if no headings the expected number of rows)

        if first_row_is_header == True:  # if selected the dictionary's keys will be created from the column titles

            for header in headings:  # iterate over the list 'headings'
                csv_dict[header.strip()] = []  # create keys in the dictionary each associated with an empty list

            """
            the strip() function will remove any surrounding whitespace characters - just makes everything a little prettier.
            """

            """
            There are two ways to do a for loop of a list:
            1: for item in list: 
                item - can be any variable - the item/variable in that iteration of the list will be contained in that varaible you defined
                
            2: for i in range(len(list)):
                i - is an integer that you can define (traditionally i, j, k, l, m etc are used but it can be anything!!!)
                
                This loop will iterate over values of i of 0 to the value put inside range() (non inclusive - this can be a little confusing -
                    i.e. for i in range(6) will iterate i from 0, 1, 2, 3, 4, 5.)
                This works really well when iterating over lists as the lists are zero-indexed!! 
                
                In-order to get the values out of the list you need to use i as an index -
                e.g. this will print all the values in a list to the console
                    for i in range(len(list)):
                        print(list[i])
            
            neither is better than the other, but the indexing option allows you to manipulate the index(i.e. look at the next value), or use it with 
            other lists of the same length - great for datasets!!!
            """

            for j in range(1, len(lines)):  # loop through each line SKIPPING THE HEADER LINE
                """
                the range function's default start point is 0. this can be changed by inputting two variables
                i.e. range(start,end) - remember range is UPTO the end value - not inclusive
                """

                line = lines[j].split(",")  # split each line into a list of strings
                for i in range(len(headings)):  # loop through each line - NOTE we're using the length of headings list

                    """
                    some rows/columns may be missing some data 
                    so if the length of a row is less than the length of the header row, we would get horrible 
                    errors!!!
                    """

                    if i < len(line):
                        """
                        you can call a value of a dictionary key like this:
                            dict_name[key]
                        
                        if the key is a string remember the quotes i.e. "key"
                        
                        you can also change the value/list/variable there too:
                        dict_name[key] = new_value
                        
                        or you can use functions on any object held in the dictonary under that key
                        i.e. if the key leads to a list - list have an append function which adds a new element to it
                            dict_name[key].append(new_value)
                        
                        """
                        csv_dict[headings[i].strip()].append(line[
                                                                 i].strip())  # add the value to the list! again using Srip() to remove whitespace characters
                    else:  # finds a missing value and replaces it with NA
                        csv_dict[headings[i].strip()].append("NA")


        else:
            """
            this is mostly the same as above, but creates keys that have numbered columns
                i.e. col 0, col 1 etc
            """
            for i in range(len(headings)):
                csv_dict["col {0}".format(i)] = []

            for j in range(len(lines)):  # don't skip the first line this time
                line = lines[j].split(",")
                for i in range(len(headings)):
                    if i < len(line):
                        csv_dict["col {0}".format(i)].append(line[i].strip())
                    else:
                        csv_dict["col {0}".format(i)].append("NA")

    return csv_dict  # this returns the dictionary as the output of the function


def csv_to_list_of_lists(file_location, first_row_is_header=True):  # creates the funtion

    with open(file_location,
              "r") as file:  # opens the file and creates a file object (don't worry if you don't understand objects yet)
        data = file.read()  # reads the whole file into memory as a single long string

        List_data = [] # define/create our output list

        # lets split the string by each line in the file
        lines = data.split(
            "\n")  # sometimes need \n\r  ! varies from csv to csv. How can we account for this? (hint: regex)

        if first_row_is_header == True:  # if selected skip first line
            for i in range(1,len(lines)):
                line = lines[i].split(",")
                List_data.append(line)
        else:
            for i in range(len(lines)):
                line = lines[i].split(",")
                List_data.append(line)

    return List_data



if __name__ == "__main__":  # don't worry about this - all this means is that the code inside this won't run unless this page is being run (which is kinda confusing right now - so don't worry)

    """
    how to use the function we created!!!
    you can use them very much like the inbuilt ones in python. function_name(input_variables)  - FYI input variables are usually called arguments
    the arguments usually have to be in the right order (in this case file_location followed by the value of first_row_is_header) but
    you can define them manually -i'll go into detail about this later
    
    
    if your function has a return value you need to give the output of it a variable name 
    
    when you call a function it will run the code inside it and if it has an return, it'll give an output
    
    
    FYI - variables/objects inside functions are not (usually) accessible outside the function 
    """
    dict_output = csv_to_dict("C:\Data\\test.csv", True)  # dict_output is the dictionary we just made!

    """
    now to test this dictionary 
    """
    for key in dict_output:
        print(key)
        print(dict_output[key])


    list_method = csv_to_list_of_lists("C:\Data\\test.csv", True)

    for i in range(len(list_method)):
        print(list_method[i])


    print(list_method[0][2])
