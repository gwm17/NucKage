#!/usr/bin/env python3

import numpy as np
import requests

class MassTable:
	def __init__(self):
		file = open("./etc/mass.txt","r")
		self.mtable = {}
		u2mev = 931.4940954
		me = 0.000548579909 #MeV
		self.etable = {}

		file.readline()
		file.readline()

		for line in file:
			entries = line.split()
			n = entries[0]
			z = entries[1]
			a = entries[2]
			element = entries[3]
			massBig = float(entries[4])
			massSmall = float(entries[5])

			key = '('+z+','+a+')'
			el_key = z
			value = (massBig+massSmall*1e-6)*u2mev - float(z)*me
			self.mtable[key] = value
			self.etable[el_key] = element
		file.close()

	def GetMass(self, z, a):
		key = '('+str(z)+','+str(a)+')'
		if key in self.mtable:
			return self.mtable[key]
		else:
			return 0

	def GetNuclearSymbol(self, z, a):
		key = str(z)
		if key in self.etable:
			return "<sup>"+str(a)+"</sup>"+self.etable[key]
		else:
			return 'none'

	def GetElementSymbol(self, z):
		key = str(z)
		if key in self.etable:
			return self.etable[key]
		else:
			return 'none'

Masses = MassTable()