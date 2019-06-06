#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__    = "Guannan Wang"
__email__     = "wangguannan2014@gmail.com"
__version__   = "1"

def isfloat(Num):
	try:
		float(Num)
		return True
	except ValueError:
		return False

def isint(Num):
	try:
		float(Num)
		int(Num)
		return True
	except ValueError:
		return False
		
def intersect(listA, listB):
    return list(set(listA) & set(listB))


def rowtolist(fin, colnum, sep, header):
	"""
	Get the a list of elements from specified column
	"""
	rowtolist = []
	fin_rowtolist = open(fin, "rU")
	for line in fin_rowtolist:
		if line not in ["\n", "\r\n"]:
			element_row = line.split(sep)[colnum - 1]
			rowtolist.append(element_row)
	if header == "Y":
		return rowtolist[1:]
	elif header == "N":
		return rowtolist

def unielement(fin, colnum, sep, header):
	"""
	Get the unique elements from specified column
	"""
	unique_element = []
	fin_unielement = open(fin, "rU")
	for line in fin_unielement:
		if line not in ["\n", "\r\n"]:
			element_row = line.split(sep)[colnum - 1]
			if element_row not in unique_element:
				unique_element.append(element_row)
	if header == "Y":
		return unique_element[1:]
	elif header == "N":
		return unique_element	

def average(values):
	"""
	Calculate the average of input vaules
	values	Values to be used for average calculation, should be a list
	"""
	numofvalues = 0 
	sumofvalues= 0.0
	for number in values:
		sumofvalues += float(number)
		numofvalues += 1
	return sumofvalues / numofvalues


def colalpha(hexcode, alpha):
	"""
	Get the hexcode of the current color with a certain transprancy, similar to "alpha" in R
	"""
	rgbcode = tuple(int(hexcode.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
	rgbcode = tuple(255 - alpha * (255 - primary) for primary in rgbcode)
	hexcode = '#%02x%02x%02x' % rgbcode
	return hexcode