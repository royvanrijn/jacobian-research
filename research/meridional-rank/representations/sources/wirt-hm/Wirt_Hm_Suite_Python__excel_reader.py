'''
    File name: excel_reader.py
    Dependencies: d4_gen.py, d5_gen.py, calc_wirt.py, sym_hm.py, d_hm.py, gauss_processor, D4_gens_pruned.txt, D5_gens_pruned.txt
    Author: Nathaniel Morrison
    Date created: 09/01/2019
    Date last modified: 03/18/2020
    Python Version: 3.7.2, but should work in any version since 2.7
'''
import xlsxwriter, sys, xlrd, gauss_processor, sym_hm, d_hm, calc_wirt
import sys
from itertools import combinations, permutations
from datetime import datetime
#In Python 2, which a lot of people still use because it's integrated into Linux and Sage and a host of other OSes and programs, the functionality provided by input() in
#Python 3 is instead provided by raw_input(). To make this program compatible with either version, define raw_input() to be input() if we're in Python 3, and use raw_input()
#in place of input() throughout the code.
if sys.version_info[0]==3:
    def raw_input(prompt):
        return input(prompt)

#Function to pull name and gauss code from Excel sheet
def knot_processor(cur_row, knot_workbook, knot_table, name_col, gauss_col):
    #Call appropriate cell and convert to string object from cell object
    name = str(knot_table.cell(cur_row, name_col))
    #Get rid of prefix
    name = name.lstrip("text:\'")
    #Get rid of suffix
    name=name.rstrip("\'")
    #Repeat for gauss code
    raw_gauss_code= str(knot_table.cell(cur_row, gauss_col))
    raw_gauss_code = raw_gauss_code.lstrip("text:\'")
    raw_gauss_code=raw_gauss_code.rstrip("\'")
    #return the name of the knot and the isolated gauss code
    return name, raw_gauss_code

#Function to create the output file
def excel_creator(excel_name, choice):
    #Create an Excel workbook object to store excel version of output
    workbook=xlsxwriter.Workbook(excel_name)
    #Add worksheet attribute within the Excel workbook object to put the data in
    worksheet = workbook.add_worksheet()
    #Enter column titles
    worksheet.write(0, 0, 'Knot Name')
    worksheet.write(0, 1, 'Gauss Notation')
    worksheet.write(0, 2, 'Seed Strand Set')
    worksheet.write(0, 3, 'Wirtinger Number')
    if choice == 2 or choice ==3:
        worksheet.write(0, 4, 'Seeed Strand Maps to Generating Set')
    #Return worksheet
    return worksheet, workbook

#Function to write data to the output file
def excel_writer(name, seed_strand_with_gauss, knot_dict, raw_gauss_code, wirt_num, sym_gen_set, Knot_Number, worksheet, workbook, choice):
    #Write data to worksheet, with the row determined by the Knot_Number accumulator
    worksheet.write(Knot_Number, 0, str(name))
    worksheet.write(Knot_Number, 1, '{' + str(raw_gauss_code) + '}')
    worksheet.write(Knot_Number, 2, str(seed_strand_with_gauss))
    worksheet.write(Knot_Number, 3, wirt_num)
    if choice == 2 or choice == 3:
        worksheet.write(Knot_Number, 4, str(sym_gen_set))
    return

def excel_main(input_name, choice):
    #Get info about input file format. Excel index start at 1 while xlrd index starts at 0, so subtract 1 from all user inputs
    name_col=int(raw_input("Please enter the number of the column containing knot names in the input file (e.g. A=1, B=2, etc): "))-1
    gauss_col=int(raw_input("Please enter the number of the column containing the gauss code in the input file: "))-1
    start_row=int(raw_input("Please enter the number of the row containing the first knot's information: "))-1
    #Get name of output and input files from user
    excel_name=raw_input("Enter desired name of excel output file (include directory and extension): ")
    #Create accumulator to let us keep track of how many knots have had a homomorphism found
    Knot_Number = 1
    #Create a workbook object and worksheet attribute
    worksheet, workbook = excel_creator(excel_name, choice)
    #Load data file.
    knot_workbook=xlrd.open_workbook(input_name)
    #Call first (and only) worksheet in workbook
    knot_table=knot_workbook.sheets()[0]
    #Get total number of rows; since Python index starts at 0, and Excel's starts at 1, subtract 1 from the value .nrows gives
    tot_rows=knot_table.nrows - 1
    #Set current row index
    cur_row = start_row
    if choice == 2:
        all_perms=sym_hm.sym_group_crafter()
    elif choice ==3:
        user_wirt=raw_input("Enter the wirtinger number of the knots you wish to analyze (only 4 and 5 are supported): ")
        while user_wirt !="4" and user_wirt!="5":
            user_wirt=raw_input("You have not entered a valid number. Please enter either 4 or 5: ")
        user_wirt=int(user_wirt)
        all_perms=d_hm.coxeter_gen_crafter(user_wirt)
    #Until the last row is reached...
    while cur_row <= tot_rows:
        name, raw_gauss_code = knot_processor(cur_row, knot_workbook, knot_table, name_col, gauss_col) #processes current row
        gauss_code = gauss_processor.process_gauss_code(raw_gauss_code) #Turns raw_gauss_code, a string, into an integer list
        knot_dict, seed_strand_set, wirt_num = calc_wirt.wirt_main(gauss_code)
        if choice == 2:
            hmorph, sym_gen_set = sym_hm.homomorphism_finder(seed_strand_set, knot_dict, wirt_num, all_perms)
            #if a homomorphism was found
            if hmorph:
                #Set up output version of the seed strand set, complete with the gauss code string of each strand
                seed_strand_with_gauss=dict()
                for strand in seed_strand_set:
                    seed_strand_with_gauss[strand]= knot_dict[strand][0]
                #write to Excel file
                excel_writer(name, seed_strand_with_gauss, knot_dict, raw_gauss_code, wirt_num, sym_gen_set, Knot_Number, worksheet, workbook, choice)
                #Iterate the accumulator
                Knot_Number = Knot_Number + 1
        elif choice == 3:
            #Set defult value of 'We have found a homomorphism' to 'no'
            hmorph = False
            #Filter out all knots with wrong wirtinger number
            if wirt_num == user_wirt:
                hmorph, cox_gen_set = d_hm.homomorphism_finder(seed_strand_set, knot_dict, wirt_num, all_perms) #Attempts to find a homomorphism
            #if a homomorphism was found
            if hmorph:
                #Set up output version of the seed strand set, complete with the gauss code string of each strand
                seed_strand_with_gauss=dict()
                for strand in seed_strand_set:
                    seed_strand_with_gauss[strand]= knot_dict[strand][0]
                #numpy arrays are great but look yucky when printed, so make 'em pretty
                cox_gen_set_writable=d_hm.cox_writer(cox_gen_set)
                #write to Excel file
                excel_writer(name, seed_strand_with_gauss, knot_dict, raw_gauss_code, wirt_num, cox_gen_set_writable, Knot_Number, worksheet, workbook, choice)
                #Iterate the accumulator
                Knot_Number = Knot_Number + 1
        else:
            #Set up output version of the seed strand set, complete with the gauss code string of each strand
            seed_strand_with_gauss=dict()
            for strand in seed_strand_set:
                seed_strand_with_gauss[strand]= knot_dict[strand][0]
            excel_writer(name, seed_strand_with_gauss, knot_dict, raw_gauss_code, wirt_num, list(), Knot_Number, worksheet, workbook, choice)
            #Iterate the accumulator
            Knot_Number = Knot_Number + 1
        cur_row+=1
    #Close the excel file
    workbook.close()
            
        
    
