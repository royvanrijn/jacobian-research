'''
    File name: sym_hm.py
    Author: Nathaniel Morrison
    Date created: 09/01/2019
    Date last modified: 03/18/2020
    Python Version: 3.7.2, but should work on any version since 2.7
'''
from itertools import combinations, permutations
from copy import copy

#Function to search for homomorphisms from the fundamental group to the symmetric group.
def homomorphism_finder(seed_strand_set, knot_dict, wirt_num,transpositions):
    #Get all permutations of the seed strand set
    perms_of_seed_strands=list(permutations(seed_strand_set))
    #Set initial value of "We have found a homomorphism" to "No"
    hmorph=False
    #Run through every generating set in the proper symmetric group
    for gen_set in transpositions[wirt_num-2]:
        #Run through every permutation of the seed strand set
        for perm in perms_of_seed_strands:
            #Map each seed strand to a generator
            sym_gen_set = generator_assign(wirt_num, perm, gen_set)
            #Map every strand to a generator
            mapping = transposition_assignment(sym_gen_set, knot_dict)
            #Determine if that mapping admits a surjective homomorphism
            hmorph = homomorphism_tester(mapping, knot_dict)
            #If it does...
            if hmorph:
                #Return the Boolean variable indicating that it does and the seed-strand-to-generator mapping that produced it
                return hmorph, sym_gen_set
    #If no hm was found across all permutations, hmorph will be false and we don't care what the program thinks the sym_gen_set is
    return hmorph, sym_gen_set

#Creates list of all generating transpositions that must be checked of each symeric group
def sym_group_crafter():
    #Defineing all the transpositions generating Sn+1 that must be checked, up to n=4
    s3a=((1,2),(1,3))
    s4a=((1,2),(1,3),(1,4))
    s4b=((1,2),(1,3),(2,4))
    s4c=((1,2),(1,3),(3,4))
    s5a=((1,2),(1,3),(1,4),(1,5))
    s5b=((1,2),(1,3),(1,4),(2,5))
    s5c=((1,2),(1,3),(1,4),(3,5))
    s5d=((1,2),(1,3),(1,4),(4,5))
    s5e=((1,2),(1,3),(2,4),(1,5))
    s5f=((1,2),(1,3),(2,4),(2,5))
    s5g=((1,2),(1,3),(2,4),(3,5))
    s5h=((1,2),(1,3),(2,4),(4,5))
    s5i=((1,2),(1,3),(3,4),(1,5))
    s5j=((1,2),(1,3),(3,4),(2,5))
    s5k=((1,2),(1,3),(3,4),(3,5))
    s5l=((1,2),(1,3),(3,4),(4,5))
    #To generate s6 generatring sets:
    #Create empty lists to house the transposition lists
    master_gen_list_s5=[]
    master_gen_list_s6=[]
    #Run through all the s5 generator sets
    for genor in (s5a,s5b,s5c,s5d,s5e,s5f,s5g,s5h,s5i,s5j,s5k,s5l):
        #Turn each tuple into a list for mutation
        gen_list_s5=list(genor)
        #Add the list version onto the big list of s5 generating sets
        master_gen_list_s5.append(gen_list_s5)
    #Run through all the s5 generating sets, list versions
    for s5_genor in master_gen_list_s5:
        #Run through the five transpositions that need to be tacked on to make s6 generating sets
        for s6_genor in [(1,6), (2,6), (3,6), (4,6), (5,6)]:
            #Create shallow copy of s5 generating set so we don't change the s5 set by accident
            #This is implemented differently in Python 2 and Python 3. So, one of these will work...
            try:
                gen_list_s6=s5_genor.copy()
            except:
                gen_list_s6=copy(s5_genor)
            #Add the current s6 transposition onto the s5 set to create an s6 generating set
            gen_list_s6.append(s6_genor)
            #Add the new s6 generating set onto the big list of s6 generating sets as a tuple
            master_gen_list_s6.append(tuple(gen_list_s6))
    #Convert the whole thing to a tuple to save memory and processing time
    s6_set=tuple(master_gen_list_s6)
    #Put all the sub-tuples into a super-tuple, because you can't have enough tuples
    transpositions=((s3a,), (s4a,s4b,s4c),(s5a,s5b,s5c,s5d,s5e,s5f,s5g,s5h,s5i,s5j,s5k,s5l), s6_set)#the extra comma in the first entry forces (s3a) to be a one-element tuple
    return transpositions


#function to assign generating transpositions to each seed strand
def generator_assign(wirt_num, perm, gen_set):
    #Defining dictionary to house generating set of transpositions
    sym_gen_set={}
    #run through each element of the permuted seed strand set
    for i in range(len(perm)):
        #Add an entry to the generating set dictionary with the letter of the seed strand from the knot dictionary as the key and the Sn+1 generator as the entry
        sym_gen_set[perm[i]]=gen_set[i]
    #return the generating set, a dictionary with the keys being the seed strand letters and the enries being the Sn+1 generators assigned to them
    return sym_gen_set

'''
Note that in the next function, I have decided to first assign a transposition to each strand, then check to make sure all crossings obey the wirtinger relations. While it
would be possible to check this while running through the assignments, I believe that for large knots this would slow the program down significantly. Large knots may require
quite a few iterations of the assignment loop, and if the program must check every already-fully-assigned crossing every iteration, the pace of calcualtion might get a bit
ridiculous. I believe it will be far faster in the long run to have two seperate loops.
'''

#Function to determine whether a given generating set gives a homomorphism to the Sn+1 group
def transposition_assignment(sym_gen_set, knot_dict):
    #Create shallow copy of the generating set dictionary
    mapping=sym_gen_set.copy()
    #This variable will indicate when all strands that can be assigned a generator have been assigned a generator
    new_assignment=True
    #while we aren't done
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
                        #Then add the second understrand to mapping, assigning it to a (possibly) new transposition
                        mapping[understrands[1]]=transpose_product(mapping[overstrand],mapping[understrands[0]])
                        #A new assignment has been made, so another iteration is necessary to make sure all possible assigments have been made
                        new_assignment=True
                    #If only the second understrand has been assigned...
                    elif understrands[1] in mapping:
                        #Then add the first understrand to mapping, assigning it to a (possibly) new transposition
                        mapping[understrands[0]]=transpose_product(mapping[overstrand],mapping[understrands[1]])
                        #A new assignment has been made, so another iteration is necessary to make sure all possible assigments have been made
                        new_assignment=True
                    #Otherwise, neither understrand has yet been assigned a transposition, so let another iteration happen
            #Otherwise, the overstrand has not been assigned a mapping yet, and there is nothing we can do yet, so let another iteration happen
    #Return mapping, a dictionary with all strands that can be assigned a transpsition as keys, and their corresponding two-tuple transpositions as entries.
    return mapping

#Function to test if the mapping given by the (aptly-named) mapping dictionary is a homoporphism from the fundamental group to the Sn+1 symmetric group
def homomorphism_tester(mapping, knot_dict):
    #If the length of mapping does not equal the length of the knot_dictionary, then something has gone badly wrong and we don't even have a valid generator assignment for
    #for each strand
    if len(mapping) != len(knot_dict):
        #Return a False value
        return False
    #Run through each entry (overstrand) in knot_dict
    for overstrand in knot_dict:
        #Scan through all its understrand pairs 
        for understrands in knot_dict[overstrand][1]:
            #If ever one of the understrands' generators (it dosen't matter which understrand, since we're working with two-cycles. Otherwise, it would) does not equal the 
            #Transpose product of the overstrand and other understrand, the Wirtinger relation does not hold for this crossing, and this is not a homomorphism.
            if set(mapping[understrands[0]]) != set(transpose_product(mapping[overstrand], mapping[understrands[1]])):
                #Return a False value
                return False
    #If the function has not been exited yet, then we have a valid homomorphism
    return True
                              
#function to find product of two transpositions, based on the relations of the Wirtinger presentation
def transpose_product(overstrands, understrands):
    #Convert each transposition tuple to a set to make computation easier
    ostrand=set(overstrands)
    ustrand=set(understrands)
    #If the two transpositions are the same...
    if ostrand == ustrand:
        #Then the product is itself
        product=overstrands
    #If the two transpositions have one element in common...
    elif len(ostrand.intersection(ustrand))== 1:
        #Then the product is the symmetric difference (the two elements that are not in both origional transpositions)
        product=tuple(ostrand.symmetric_difference(ustrand)) #note that the product is converted back to a tuple so that it can be iterated over later
    #otherwise, the two transpositions are entirely disjoint
    else:
        #In this case, the product is the understrand
        product=understrands
    #Return the product, a two-tuple
    return product
