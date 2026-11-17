'''
    File name: d_hm.py
    Dependencies: d4_gen.py, d5_gen.py, D4_gens_pruned.txt, D5_gens_pruned.txt
    Author: Nathaniel Morrison
    Date created: 09/01/2019
    Date last modified: 11/03/2019
    Python Version: 3.7.2
'''
from itertools import combinations, permutations
import numpy as np
import d4_gen
import d5_gen
from copy import copy

#Function to search for homomorphisms from the fundamental group to the coxeter group.
def homomorphism_finder(seed_strand_set, knot_dict, wirt_num,cox_perms):
    #n will give which generating set to use in this iteration
    n=0
    #Set initial value of "We have found a homomorphism" to "No"
    hmorph=False
    #While we still have permutations left to check and none have worked so far
    #n_set=set()
    while n<len(cox_perms) and hmorph == False:
        #Map seed strands to generators
        cox_gen_set=generator_assign(seed_strand_set, wirt_num, cox_perms[n])
        #Map all the other strands to coexter elements as per the Wirtinger relations
        mapping=transposition_assignment(cox_gen_set, knot_dict)
        #Check whether the final mapping is a surjective homomorphism by checking the Wirtinger realtion at each crossing
        hmorph=homomorphism_tester(mapping, knot_dict)
        n+=1
    #Return the Boolean variable stating whether there is a valid homomorphism or not, and the generating set
    return hmorph, cox_gen_set

#Function to create a list of all generators of the requested D group that must be checked
def coxeter_gen_crafter(user_wirt):
    if user_wirt ==4:
        #Instantiate object
        d4gen_inst=d4_gen.DGen()
        #Call main() method within d4_gen module
        D4_gen_sets=d4gen_inst.dgen_main()
        #Create a list to house the permuted generating sets
        D4_gen_perms_lis=list()
        #Let the uesr know it's working
        print("Computing permutations. This may take some time...")
        #Run through all the generating sets recieved from the module
        for gen_set in D4_gen_sets:
            #Find all permutations of each set and add them to the list of all permuted sets that must be checked
            D4_gen_perms_lis=D4_gen_perms_lis+list(permutations(gen_set))
        #Convert list back to tuple because we're done mutating it
        D4_gen_perms=tuple(D4_gen_perms_lis)
        #Let the user know it's worked
        print("Permutations computed!")
        #Assign the list to the output variable.
        all_perms=D4_gen_perms
    #Do same if wirtinger number is 5
    elif user_wirt ==5:
        d5gen_inst=d5_gen.D5Gen()
        D5_gen_sets=d5gen_inst.d5gen_main()
        D5_gen_perms_lis=list()
        print("Computing permutations. This may take some time...")
        for gen_set in D5_gen_sets:
            D5_gen_perms_lis=D5_gen_perms_lis+list(permutations(gen_set))
        D5_gen_perms=tuple(D5_gen_perms_lis)
        print("Permutations computed!")
        all_perms=D5_gen_perms
    return all_perms


#Function to assign generators to each seed strand
def generator_assign(seed_strand_set, wirt_num, cox_perm):
    #Defining dictionary to house generating sets
    cox_gen_set={}
    #Run through each element of this iteration's coxeter generators
    for i in range(len(cox_perm)):
        #Add an entry to the generating set dictionary with the letter of the seed strand from the knot dictionary as the key and the generator as the entry
        cox_gen_set[seed_strand_set[i]]= cox_perm[i]
    #Return the generating set, a dictionary with the keys being the seed strand letters and the enries being the generators assigned to them
    return cox_gen_set

'''
Note that in the next function, I have decided to first assign a generator to each strand, then check to make sure all crossings obey the wirtinger relations. While it
would be possible to check this while running through the assignments, I believe that for large knots this would slow the program down significantly. Large knots may require
quite a few iterations of the assignment loop, and if the program must check every already-fully-assigned crossing every iteration, the pace of calcualtion might get a bit
ridiculous. I believe it will be far faster in the long run to have two seperate loops.
'''

#Function to determine whether a given generating set gives a homomorphism to the coxeter group. It's called "transposition" assignment because that's what I called it when I
#was writing the symetric group code and I couldn't be bothered to track down everywhere the function name appears and change it.
def transposition_assignment(cox_gen_set, knot_dict):
    #Create shallow copy of the generating set dictionary
    mapping=cox_gen_set.copy()
    #This variable will indicate when all strands that can be assigned a generator have been assigned a generator
    new_assignment=True
    #While we aren't done
    while new_assignment:
        #This iteration will be the last one unless told otherwise
        new_assignment=False
        #Scan through each overstrand in the knot dictionary
        for overstrand in knot_dict:
            #If the overstrand has been assigned a mapping...
            if overstrand in mapping:
                #Scan through all its understrand pairs (there may be more than one if the strand is the overstrand of more than one crossing)
                for understrands in knot_dict[overstrand][1]:
                    #If both understrands of this particular crossing have been assigned...
                    if understrands[0] in mapping and understrands[1] in mapping:
                        #This crossing has been fully assigned and we're done here
                        pass
                    #If only the first understrand has been assigned...
                    elif understrands[0] in mapping:
                        #Then add the second understrand to mapping, assigning it to a (possibly) new group element
                        mapping[understrands[1]]=transpose_product(mapping[overstrand],mapping[understrands[0]])
                        #A new assignment has been made, so another iteration is necessary to make sure all possible assigments have been made
                        new_assignment=True
                    #If only the second understrand has been assigned...
                    elif understrands[1] in mapping:
                        #Then add the first understrand to mapping, assigning it to a (possibly) new group element
                        mapping[understrands[0]]=transpose_product(mapping[overstrand],mapping[understrands[1]])
                        #A new assignment has been made, so another iteration is necessary to make sure all possible assigments have been made
                        new_assignment=True
                    #Otherwise, neither understrand has yet been assigned a group element, so let another iteration happen
            #Otherwise, the overstrand has not been assigned a mapping yet, and there is nothing we can do yet, so let another iteration happen
    #Return mapping, a dictionary with all strands that can be assigned a group element as keys, and their corresponding numpy array Coxeter group elements as entries.
    return mapping

#Function to test if the mapping given by the (aptly-named) mapping dictionary is a homoporphism from the fundamental group to the Coxeter group
def homomorphism_tester(mapping, knot_dict):
    #If the length of mapping does not equal the length of the knot_dictionary, then something has gone badly wrong and we don't even have a valid generator assignment for
    #for each strand
    if len(mapping) != len(knot_dict):
        #Return a False value
        return False
    #Run through each entry (overstrand) in knot_dict
    for overstrand in knot_dict:
        #scan through all its understrand pairs 
        for understrands in knot_dict[overstrand][1]:
            #If ever one of the understrands' generators (it dosen't matter which understrand, since all the group elements are conveniently Hermitian. Otherwise, it would)
            #does not equal the transpose product of the overstrand and other understrand, the Wirtinger relation does not hold for this crossing, and this is not a
            #homomorphism.
            if not np.array_equal(mapping[understrands[0]], transpose_product(mapping[overstrand], mapping[understrands[1]])):
                #return a False value
                return False
    #If the function has not been exited yet, then we have a valid homomorphism
    return True
                              
#Function to find product of three matrices as per the Wirtinger relations. It's called "transpose" product because that's what it was called when I was doing the symmetric
#groups and I couldn't be bothered to track down and change every spot where the name appears. Nothing is actually being transposed here.
def transpose_product(ostrand, ustrand):
    #Just use numpy methods.
    semiproduct=np.dot(ustrand,ostrand)
    product=np.dot(ostrand,semiproduct)
    return product

#Function to convert array objects into better-looking, formatted string objects for writing to Excel
def cox_writer(cox_gen_set):
    #Create a new dictionary of beutified matrices
    cox_gen_set_writable=dict()
    #Gather a list of the dictionary keys
    cox_keys=list(cox_gen_set.keys())
    #Run through each key
    for key in cox_keys:
        #Create an empty string to house the new matrix
        matrix=''
        #Run through all the rows of the numpy array, which are themselves arrays, after converting the larger array into a tuple
        for array_row in tuple(cox_gen_set[key]):
            #Make it look pretty. Or at least less ugly
            array_row_str=str(tuple(array_row))
            matrix=matrix+array_row_str+"|"
        matrix=matrix.rstrip("|")
        #Reload The Matrix (ba-dum tish!) into the dictionary
        cox_gen_set_writable[key]=matrix
    return cox_gen_set_writable
