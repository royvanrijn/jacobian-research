'''
    Module name: calc_wirt.py
    Author: Paul Villanueva
    Edited By: Nathaniel Morrison
    Date created: 12/8/2018
    Date last modified: 03/18/2020
    Python Version: 3.7.2, but should work with any verison since 2.7
'''

from itertools import combinations, permutations
from copy import copy

#Function to manage the creation of the knot dictionary
def create_knot_dictionary(gauss_code):   
    #Isolate each strand of the knot diagram from the gauss code
    strands_dict = find_strands(gauss_code)
    #Find the crossing information of each strand
    knot_dict = find_crossings(strands_dict, gauss_code)
    #Return the knot dictionary, a dictionary in which each key is a strand label and each entry is a 2-element list, the first element being the subset of the gauss
    #code defining the strand and the second being a list of the crossings for which the strand is the overstrand
    return knot_dict

#Function to find and label each strand in the knot diagram from the gauss code. Each strand is a subinterval of the gauss code lying between two negative numbers
def find_strands(gauss_code):
    #Creates a set to store the strands
    strand_set = set()
    i = 0
    #Do this loop indefinitely
    while True:
        #If the ith entry of the gauss_code list is negative...
        if gauss_code[i] < 0:
            #Declare this to be the beginning of the strand
            beginning = i
            #Iterate i mod len(gauss_code); that way, if i > len(gauss_code), the program "wraps around" back to the first entry
            i = (i + 1) % len(gauss_code)
            #Keep running though the gauss code as long as the numbers are positive
            while gauss_code[i] > 0:
                #If the current entry is positive, move on to the next one
                i = (i + 1) % len(gauss_code)
            #If the beginning entry is later in the list than the next negative one (e.g. the program has wrapped around back to the start of the list)...
            if beginning > i:
                #Then the strand is all values from the beginning value to the end of the list,...
                new_strand = gauss_code[beginning:]
                #... and all values from the first entry to the first negative value
                for k in range(i + 1):
                    new_strand.append(gauss_code[k])
                #Convert the new_strand list to a tuple; takes up less memory
                new_strand = tuple(new_strand)
            #If the program has not had to wrap around to the start of the list...
            else:
                #Then the new_strand is just the gauss code from the beginning value to the negative value
                new_strand = tuple(gauss_code[beginning:i+1])
            #If this strand has not previously been found...
            if new_strand not in strand_set:
                #Then add the strand to the set of all strands
                strand_set.add(new_strand)
            #If it has already been found, exit the while loop
            else:
                break
        #If the ith entry is not negative, look at the next one
        else:
            i = (i + 1) % len(gauss_code)
    #Define a list of all capital letters of the English alphabet
    letter_list = list(map(chr, range(65, 91)))
    #Expanding strand labeling to work for knots up to 702 strands by including labels of form 'AA', 'AB', etc, past A-Z:
    #If more strand labels are needed(technically it's less than one strand per two crossings, but close enough)...
    if (len(gauss_code))/2 >26:
        #Determine how many more strands are needed (again, technically much less than this, but close enough. This piece of code dosen't take that long to iterate)
        diff = len(gauss_code) - 26
        #Determine how many more sets of 26 strand labels are needed
        needed_labels = diff/26
        j= 65
        #While we still don't have enough sets of 26...
        while j-65 < needed_labels:
            l = 65
            #While l is less than 90 (which is the ascii code for Z)
            while l <= 90:
                #Append the appropiate doble letter 'j'+'l' to the list of labels
                letter_list.append(chr(j) + chr(l))
                l=l+1
            j=j+1
    #Create a dictionary to house the strands
    strands_dict = dict()
    #Enumerate retuns a list of two-tuples where the first entry is an integer i running from 1 to len(strand_set) and the second is the ith entry
    for i, strand in enumerate(strand_set):
        #Assign a letter label to each strand in strand_set
        strands_dict[letter_list[i]] = [strand, []]
    #Return the strand dictionary
    return strands_dict

#Function to finish the knot dictoionary by adding crossing information
#Runs through each strand in strand_dict. At each one, runs through each positive number in its gauss code. Looks for the two strands that begin or end with the negative
#version of that number. These are the understands of that crossing. Once these have been found, add them to the list of crossings associated with that overstrand, and
#repeat until all understrands of all crossings have been found.
def find_crossings(knot_dict, gauss_code):
    #Run through each key (strand) in knot_dict
    for key_outer in knot_dict:
        #Run through each element of the gauss code corresponding to the strand
        for under in knot_dict[key_outer][0]:
            #If the element is not negatve, this strand may be an overstrand of a new crossing, so declare found1, found2 false to move onto next step
            if under > 0:
                found1, found2 = False, False
                #Run through each key (strand) in knot_dict
                for key_inner in knot_dict:
                    #If both the understands are in the knot dictionary already, exit loop
                    if found1 and found2:
                        break
                    #If at least one understrand isn't...
                    else:
                        #If the the current gauss code element is minus the first element in the gauss code of one of the entries in the dictionary...
                        if knot_dict[key_inner][0][0] == -under:
                            #Store this value as one of the understrands
                            under1 = key_inner
                            #Say that one of two understrands has been found
                            found1 = True
                        #If the the current gauss code element is minus the last element in the gauss code of one of the entries in the dictionary...
                        if knot_dict[key_inner][0][-1] == -under:
                            #Store this value as one of the understrands
                            under2 = key_inner
                            #Say that one of two understrands has been found
                            found2 = True
                #Add this crossing the the list associated with the overstrand
                knot_dict[key_outer][1].append((under1, under2))
    #Return the completed dictionary
    return knot_dict

'''
A note on what the knot dictionary is:
The output knot dictionary is of the form
            d_k = {
                s_i: [(gauss_subseq), [c_1, c_2,...]]
                .
                .
                .
            },
where s_i is the name of the strand, gauss_subseq is a tuple representing the subsequence of the Gauss code
corresponding to the the strand, and the c_i represent the crossings that s_n are over.  The c_i are tuples 
(s_i_1, s_i_2), where s_i_1 and s_i_2 are the names of the strands that s_i is over.      
'''

#Function to perform coloring moves as per the definition of the wirtinger number using generated seed strands to determine if they produce a valid, full coloring.
def is_valid_coloring(seed_strands, knot_dict):
    #Turn the seed stands into a set to weed out repeats
    seed_strands = set(seed_strands)
    #Shallow copy the set, so that alterations to seed_strands will not be made to colored_set (and vice versa)
    #Shallow copies are implemented differently in Python 2 and Python 3. So, one of these two will work...
    try:
        colored_set = seed_strands.copy()
    except:
        colored_set = copy(seed_strands)
    #Initialize new coloring flag
    new_coloring = True
    #While a new coloring was made during the last iteration...
    while new_coloring:
        #Assume no new coloring could be done
        new_coloring = False
        try:
            colored_set_copy = colored_set.copy()
        except:
            colored_set_copy = copy(colored_set)
        #Run through each strand in a shallow copy of colored_set. These are all the strands that are currently colored
        for strand in colored_set_copy:
            #Run through each crossing two-tuple associated with the given strand
            for crossing in knot_dict[strand][1]:
                #If one of the understrands has not been colored...
                if crossing[0] not in colored_set or crossing[1] not in colored_set:
                    #If one of the understrands has been colored...
                    if crossing[0] in colored_set or crossing[1] in colored_set:
                        #Add the two understrands to the set of colored strands (since it is a set, only the one that was not previously colored will be added)
                        colored_set.update(crossing)
                        #Say that a new coloring has been found, so the loop may repeat
                        new_coloring = True
    #Once all strands that can be colored have been colored, test if all the strands in the knot have been colored
    if colored_set == set(knot_dict.keys()):
        #If they have, return true
        return True
    #If not, return false
    return False

#Function generates seed strands and calls is_valid_coloring to test if they give a full coloring. Finds the smallest set of seed strands and thus w(D) for the given diagram
def calc_wirt_info(knot_dict):
    #Need at least two seed strands (for anything but the unknot)
    n = 2
    #Can have at most all but one strands as seed strands
    while n < len(knot_dict):
        #Run through every possible set (combination, not permutation) of n strands of the knot to use as seed strands
        for seed_strands in combinations(knot_dict, n):
            #If that set leads to the whole knot being colored...
            if is_valid_coloring(seed_strands, knot_dict):
                #Return this set of seed strands. Since n is starting from the minimum (2), this is a minimal set, and its length is w(D)
                return (seed_strands, n)
        #If no n-length subsets provide a good generating set, iterate and move on to all n+1 length sets
        n += 1
    # If control passes the above while loop, then a valid coloring with less than n - 1 strands
    # was not found.  In that case, return n - 1 colorability and an arbitrary subset
    # of n - 1 strands.
    return (set(knot_dict.keys()).pop, n)

def wirt_main(gauss_code):
    knot_dict = create_knot_dictionary(gauss_code) #Generates a knot dictionary from gauss code
    seed_strand_set, wirt_num = calc_wirt_info(knot_dict) #Determines the smallest length seed strand set and its length, which is w(d)
    return knot_dict, seed_strand_set, wirt_num

