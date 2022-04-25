#!/usr/bin/env python

import sys
from qtpy.QtWidgets import QApplication, QWidget, QMainWindow
from qtpy.QtWidgets import QLabel, QMenuBar, QAction
from qtpy.QtWidgets import QHBoxLayout, QVBoxLayout, QGridLayout, QGroupBox
from qtpy.QtWidgets import QPushButton, QButtonGroup, QRadioButton
from qtpy.QtWidgets import QSpinBox, QDoubleSpinBox, QComboBox
from qtpy.QtWidgets import QDialog, QFileDialog, QDialogButtonBox
from qtpy.QtWidgets import QTableWidget, QTableWidgetItem
from qtpy.QtWidgets import QLineEdit, QTabWidget, QFormLayout, QListWidget, QListWidgetItem
from qtpy import QtWidgets
from qtpy import QtGui
from qtpy.QtCore import Signal
from NucData import Masses
import enum

class Target:
	def __init__(self, z=[], s=[], thickness=0.0):
		self.Z = z
		self.S = s
		self.thickness = thickness

	def __str__(self):
		val = ""
		for i in range(len(self.Z)):
			val += Masses.GetElementSymbol(self.Z[i]) + "<sub>"+str(self.S[i])+"</sub>"
		return val

class Nucleus:
	def __init__(self, z=0, a=0):
		self.Z = z
		self.A = a
		self.symbol = Masses.GetNuclearSymbol(self.Z, self.A)
		self.mass = Masses.GetMass(self.Z, self.A)

	def SetNucleus(self, z, a):
		self.Z = z
		self.A = a
		self.symbol = Masses.GetNuclearSymbol(self.Z, self.A)

	def __add__(self, other):
		return Nucleus(self.Z + other.Z, self.A + other.A)

	def __sub__(self, other):
		return Nucleus(z=(self.Z - other.Z), a=(self.A - other.A))

	def __eq__(self, other):
		if self.Z == other.Z and self.A == other.A:
			return True
		else:
			return False

	def __str__(self):
		return self.symbol

class ReactorType(enum.Enum):
	Invalid = 0
	Decay = 1
	Reaction = 2

class Reactor:
	def __init__(self, nuclei=[], residualExMean=0.0, residualExSigma=0.0, beamEnergyMean=0.0, beamEnergySigma=0.0):
		self.SetNuclei(nuclei)
		self.residualExMean = residualExMean
		self.residualExSigma = residualExSigma
		self.beamEnergyMean = beamEnergyMean
		self.beamEnergySigma = beamEnergySigma

	def SetNuclei(self, nuclei):
		self.nuclei = nuclei
		if len(nuclei) == 2:
			self.nuclei.append(self.nuclei[0] - self.nuclei[1])
		elif len(nuclei) == 3:
			self.nuclei.append(self.nuclei[0] + self.nuclei[1] - self.nuclei[2])
		else:
			self.nuclei = []
		self.SetType()

	def SetType(self):
		if len(self.nuclei) == 4:
			self.type = ReactorType.Reaction
		elif len(self.nuclei) == 3:
			self.type = ReactorType.Decay
		else:
			self.type = ReactorType.Invalid

	def CheckThreshold(self, targetEx=0.0):
		if self.type is ReactorType.Reaction:
			Q = targetEx + self.nuclei[0].mass + self.nuclei[1].mass - self.nuclei[2].mass - self.nuclei[3].mass - self.residualExMean
			Ethresh = -Q * (self.nuclei[2].mass + self.nuclei[3].mass) / (self.nuclei[2].mass + self.nuclei[3].mass - self.nuclei[1].mass)
			if self.beamEnergyMean < Ethresh:
				return False
			else:
				return True
		elif self.type is ReactorType.Decay:
			Q = targetEx + self.nuclei[0].mass - self.nuclei[1].mass - self.nuclei[2].mass - self.residualExMean
			if Q < 0:
				return False
			else:
				return True
		else:
			return False

	def CheckNuclei(self):
		for nuc in self.nuclei:
			if nuc.symbol == "none":
				return False
		return True

	def GetTarget(self):
		return self.nuclei[0]

	def GetResidual(self):
		if self.type == ReactorType.Reaction:
			return self.nuclei[3]
		elif self.type == ReactorType.Decay:
			return self.nuclei[2]
		else:
			return Nucleus()

	def __str__(self):
		if self.type == ReactorType.Decay:
			return str(self.nuclei[0])+"->"+str(self.nuclei[1])+"+"+str(self.nuclei[2])
		elif self.type == ReactorType.Reaction:
			return str(self.nuclei[0])+"("+str(self.nuclei[1])+str(self.nuclei[2])+")"+str(self.nuclei[3])

	def StringNoTarg(self):
		if self.type == ReactorType.Decay:
			return "->"+str(self.nuclei[1])+str(self.nuclei[2])
		elif self.type == ReactorType.Reaction:
			return "("+str(self.nuclei[1])+str(self.nuclei[2])+")"+str(self.nuclei[3])



class ReactorChain:
	def __init__(self, reactors=[], target=None):
		self.reactors = reactors
		self.target = target
		
	def AddReactor(self, nuclei, residualExMean=0.0, residualExSigma=0.0, beamEnergyMean=0.0, beamEnergySigma=0.0):
		self.reactors.append(Reactor(nuclei, residualExMean, residualExSigma, beamEnergyMean, beamEnergySigma))

	def VerifyChain(self):
		if len(self.reactors) == 0:
			return False

		if self.reactors[0].CheckNuclei() == False or self.reactors[0].CheckThreshold() == False:
			return False

		for i in range(1, len(self.reactors)):
			if self.reactors[i].CheckNuclei() == False or not self.reactors[i].CheckThreshold(self.reactors[i-1].residualExMean) or not self.reactors[i].GetTarget() == self.reactors[i-1].GetResidual():
				return False
		return True

	def __str__(self):
		val = str(self.reactors[0])
		for i in range(1, len(self.reactors)):
			val += self.reactors[i].StringNoTarg()
		return val

class DetectorArray:
	def __init__(self, detectors=[], detector_args=[]):
		self.detectors = detectors
		self.detector_args = detector_args

	def AddDetector(self, detname, detargs):
		self.detectors.append(detname)
		self.detector_args.append(detargs)


class Config:
	def __init__(self):
		self.reactorChains = []
		self.array = DetectorArray()
		self.samples = 0
		self.outputfile = ""

	def Write(self, configfile):
		try:
			with open(configfile, "w") as file :
				file.write("begin_simulator\n")
				file.write("\t" + self.outputfile + "\n")
				file.write("\t" + str(self.samples) + "\n")
				for chain in self.reactorChains:
					file.write("\tbegin_reactorchain\n")
					for reactor in chain.reactors:
						file.write("\t\tbegin_reactor\n")
						file.write("\t\t\t" + str(reactor.residualExMean) + "\n")
						file.write("\t\t\t" + str(reactor.residualExSigma) + "\n")
						file.write("\t\t\t" + str(reactor.beamEnergyMean) + "\n")
						file.write("\t\t\t" + str(reactor.beamEnergySigma) + "\n")
						file.write("\t\t\tbegin_nuclei\n")
						for i in range(len(reactor.nuclei)-1):
							file.write("\t\t\t\t" + str(reactor.nuclei[i].Z) + " " + str(reactor.nuclei[i].A) + "\n")
						file.write("\t\t\tend_nuclei\n")
						file.write("\t\tend_reactor\n")
	
					file.write("\t\tbegin_target\n")
					file.write("\t\t\t" + str(chain.target.thickness) + "\n")
					file.write("\t\t\tbegin_elements\n")
					for i in range(len(chain.target.Z)):
						file.write("\t\t\t\t" + str(chain.target.Z[i]) + " " + str(chain.target.S[i]) + "\n")
					file.write("\t\t\tend_elements\n")
					file.write("\t\tend_target\n")
					file.write("\tend_reactorchain\n")
	
				file.write("\tbegin_detectorarray\n")
				for i in range(len(self.array.detectors)):
					file.write("\t\t" + self.array.detectors[i] + " " + self.array.detector_args[i] + "\n")
				file.write("\tend_detectorarray\n")
				file.write("end_simulator\n")
		except IOError as error:
			print("Error opening role file ", configfile)

class RoleGUI(QMainWindow):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("RoleGUI")
		self.roleConfig = Config()
		self.targets = []
		self.reactors = []
		self.chains =  []

		self.generalLayout = QVBoxLayout()
		self.centralWidget = QTabWidget(self)
		self.setCentralWidget(self.centralWidget)
		self.centralWidget.setLayout(self.generalLayout)

		self.chainTab = QWidget(self.centralWidget)
		self.chainLayout = QVBoxLayout()
		self.chainTab.setLayout(self.chainLayout)
		self.reactorTab = QWidget(self.centralWidget)
		self.reactorLayout = QVBoxLayout()
		self.reactorTab.setLayout(self.reactorLayout)
		self.targetTab = QWidget(self.centralWidget)
		self.targetLayout = QVBoxLayout()
		self.targetTab.setLayout(self.targetLayout)

		self.centralWidget.addTab(self.chainTab, "Reaction Chains")
		self.centralWidget.addTab(self.reactorTab, "Reactions")
		self.centralWidget.addTab(self.targetTab, "Targets")

		self.CreateChainTab()
		self.CreateReactorTab()
		self.CreateTargetTab()
		self.show()

	def CreateChainTab(self):
		self.roleGroup = QGroupBox("Role File", self.chainTab)
		roleLayout = QHBoxLayout()
		self.roleNameInput = QLineEdit("./roles/test.role")
		self.roleOpenButton = QPushButton("Open", self.roleGroup)
		self.roleOpenButton.clicked.connect(self.HandleRoleOpenButton)
		roleLayout.addWidget(self.roleNameInput)
		roleLayout.addWidget(self.roleOpenButton)
		self.roleGroup.setLayout(roleLayout)
		self.chainLayout.addWidget(self.roleGroup)

		self.selectedGroup = QGroupBox("Simulation Options", self.chainTab)
		selectedLayout = QHBoxLayout()
		self.selectedGroup.setLayout(selectedLayout)

		self.chainGroup = QGroupBox("Reaction Chains", self.selectedGroup)
		chainLayout = QVBoxLayout()
		self.chainTable = QTableWidget(self.chainGroup)
		self.chainTable.setColumnCount(2)
		self.chainTable.setHorizontalHeaderLabels(["Chain Equation","Target"])
		selectedLayout.addWidget(self.chainTable)
		self.chainTable.resizeColumnsToContents()
		chainLayout.addWidget(self.chainTable)
		self.chainButton = QPushButton("Add", self.chainGroup)
		self.chainButton.clicked.connect(self.HandleChainButton)
		chainLayout.addWidget(self.chainButton)
		self.chainGroup.setLayout(chainLayout)
		selectedLayout.addWidget(self.chainGroup)

		self.parametersGroup = QGroupBox("Parameters", self.selectedGroup)
		parLayout = QVBoxLayout()
		samplesLabel = QLabel("Number of Samples")
		self.samplesInput = QSpinBox(self.parametersGroup)
		self.samplesInput.setRange(0, 9000000000)
		parLayout.addWidget(samplesLabel)
		parLayout.addWidget(self.samplesInput)
		self.focalPlaneButton = QRadioButton("Focal Plane Detector", self.parametersGroup)
		self.focalPlaneButton.setAutoExclusive(False)
		fpBfieldLabel = QLabel("Focal Plane BField (kG)")
		self.fpBfieldInput = QDoubleSpinBox(self.parametersGroup)
		self.fpBfieldInput.setRange(0.0, 16.0)
		self.fpBfieldInput.setDecimals(4)
		fpAngleLabel = QLabel("Focal Plane Angle (deg)")
		self.fpAngleInput = QDoubleSpinBox(self.parametersGroup)
		self.fpAngleInput.setRange(0.0, 180.0)
		self.sabreButton = QRadioButton("SABRE Detector", self.parametersGroup)
		self.sabreButton.setAutoExclusive(False)
		parLayout.addWidget(self.focalPlaneButton)
		parLayout.addWidget(fpBfieldLabel)
		parLayout.addWidget(self.fpBfieldInput)
		parLayout.addWidget(fpAngleLabel)
		parLayout.addWidget(self.fpAngleInput)
		parLayout.addWidget(self.sabreButton)
		self.parametersGroup.setLayout(parLayout)
		selectedLayout.addWidget(self.parametersGroup)

		self.chainLayout.addWidget(self.selectedGroup)

		self.nameGroup = QGroupBox("Simulation Output File", self.chainTab)
		nameLayout = QHBoxLayout()
		self.outNameInput = QLineEdit("./test.root")
		self.outputOpenButton = QPushButton("Open", self.nameGroup)
		self.outputOpenButton.clicked.connect(self.HandleOutputOpenButton)
		nameLayout.addWidget(self.outNameInput)
		nameLayout.addWidget(self.outputOpenButton)
		self.nameGroup.setLayout(nameLayout)
		self.chainLayout.addWidget(self.nameGroup)

		self.writeButton = QPushButton("Write", self.chainTab)
		self.writeButton.clicked.connect(self.HandleWriteButton)
		self.chainLayout.addWidget(self.writeButton)

	def HandleChainButton(self):
		dia = ChainDialog(parent=self, reactors=self.reactors, targets=self.targets)
		dia.new_chain.connect(self.AddChain)
		dia.exec()
		return

	def AddChain(self, chain):
		self.chains.append(chain)
		self.UpdateChainTable()

	def UpdateChainTable(self):
		self.chainTable.setRowCount(len(self.chains))
		for i in range(len(self.chains)):
			self.chainTable.setCellWidget(i, 0, QLabel(str(self.chains[i])))
			self.chainTable.setCellWidget(i, 1, QLabel(str(self.chains[i].target)))
		self.chainTable.resizeColumnsToContents()
		self.chainTable.resizeRowsToContents()

	def HandleRoleOpenButton(self):
		file = QFileDialog.getSaveFileName(self, "Role File", "./roles", "Role File (*.role)")
		if file[0]:
			self.roleNameInput.setText(file[0])
		return

	def HandleOutputOpenButton(self):
		file = QFileDialog.getSaveFileName(self, "Output File", "./", "ROOT File (*.root)")
		if file[0]:
			self.outNameInput.setText(file[0])
		return

	def HandleWriteButton(self):
		new_config = Config()
		new_config.reactorChains = self.chains
		new_config.samples = self.samplesInput.value()
		new_config.outputfile = self.outNameInput.text()

		if self.focalPlaneButton.isChecked():
			sps_args = str(self.fpAngleInput.value())+" "+str(self.fpBfieldInput.value())
			new_config.array.AddDetector("focalplane", sps_args)
		if self.sabreButton.isChecked():
			new_config.array.AddDetector("sabre", "")

		new_config.Write(self.roleNameInput.text())
		return

	def CreateReactorTab(self):
		self.reactorGroup = QGroupBox("Reactions", self.reactorTab)
		reactorLayout = QVBoxLayout()
		self.reactorTable = QTableWidget(self.reactorGroup)
		self.reactorTable.setColumnCount(5)
		self.reactorTable.setHorizontalHeaderLabels(["Equation","Ex Mean (MeV)","Ex Sigma (MeV)","BeamE Mean (MeV)","BeamE Sigma (MeV)"])
		reactorLayout.addWidget(self.reactorTable)
		self.reactorGroup.setLayout(reactorLayout)
		self.reactorLayout.addWidget(self.reactorGroup)
		self.reactorTable.resizeColumnsToContents()

		self.reactorButton = QPushButton("Add", self.reactorGroup)
		self.reactorButton.clicked.connect(self.HandleReactorButton)
		self.reactorLayout.addWidget(self.reactorButton)
		#self.reactorTable.cellDoubleClicked.connect(self.HandleUpdateReactor)

	def HandleReactorButton(self):
		dia = ReactorDialog(self)
		dia.new_reactor.connect(self.AddReactor)
		dia.exec()
		return

	def AddReactor(self, reactor):
		self.reactors.append(reactor)
		self.UpdateReactorTable()

	def UpdateReactorTable(self):
		self.reactorTable.setRowCount(len(self.reactors))
		for i in range(len(self.reactors)):
			self.reactorTable.setCellWidget(i, 0, QLabel(str(self.reactors[i])))
			self.reactorTable.setItem(i, 1, QTableWidgetItem(str(self.reactors[i].residualExMean)))
			self.reactorTable.setItem(i, 2, QTableWidgetItem(str(self.reactors[i].residualExSigma)))
			self.reactorTable.setItem(i, 3, QTableWidgetItem(str(self.reactors[i].beamEnergyMean)))
			self.reactorTable.setItem(i, 4, QTableWidgetItem(str(self.reactors[i].beamEnergySigma)))
		self.reactorTable.resizeColumnsToContents()
		self.reactorTable.resizeRowsToContents()

	def CreateTargetTab(self):
		self.targetGroup = QGroupBox("Targets", self.targetTab)
		targetLayout = QVBoxLayout()
		self.targetTable = QTableWidget(self.targetGroup)
		self.targetTable.setColumnCount(2)
		self.targetTable.setHorizontalHeaderLabels(["Composition","Thickness (ug/cm^2)"])
		targetLayout.addWidget(self.targetTable)
		self.targetGroup.setLayout(targetLayout)
		self.targetLayout.addWidget(self.targetGroup)
		self.targetTable.resizeColumnsToContents()

		self.targButton = QPushButton("Add", self.targetGroup)
		self.targButton.clicked.connect(self.HandleTargetButton)
		self.targetLayout.addWidget(self.targButton)
		#self.targetTable.cellDoubleClicked.connect(self.HandleUpdateTarget)

	def HandleTargetButton(self):
		dia = TargetDialog(self)
		dia.new_target.connect(self.AddTarget)
		dia.exec()
		return

	def AddTarget(self, target):
		self.targets.append(target)
		self.UpdateTargetTable()

	def UpdateTargetTable(self):
		self.targetTable.setRowCount(len(self.targets))
		for i in range(len(self.targets)):
			self.targetTable.setCellWidget(i, 0, QLabel(str(self.targets[i])))
			self.targetTable.setItem(i, 1, QTableWidgetItem(str(self.targets[i].thickness)))
		self.targetTable.resizeColumnsToContents()
		self.targetTable.resizeRowsToContents()

class ChainDialog(QDialog):
	new_chain = Signal(ReactorChain)
	def __init__(self, parent=None, reactors=[], targets=[], chain=None):
		super().__init__(parent)
		self.setWindowTitle("Add A Chain")

		self.reactors = reactors
		self.targets = targets
		self.selectReactors = []
		self.selectedTarget = None

		QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		self.buttonBox = QDialogButtonBox(QBtn)
		self.buttonBox.accepted.connect(self.accept)
		self.buttonBox.accepted.connect(self.SendChain)
		self.buttonBox.rejected.connect(self.reject)
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

		self.CreateChainInputs()
		if chain is not None:
			self.SetInitialValues(chain)
		self.layout.addWidget(self.buttonBox)

	def CreateChainInputs(self):
		self.listGroup = QGroupBox("Selections", self)
		listLayout = QHBoxLayout()
		self.listGroup.setLayout(listLayout)

		self.reactorInGroup = QGroupBox("Input Reactions", self)
		reacInLayout = QVBoxLayout()
		self.reactorList = QListWidget(self.reactorInGroup)
		for reactor in self.reactors:
			item = QListWidgetItem()
			self.reactorList.addItem(item)
			self.reactorList.setItemWidget(item, QLabel(str(reactor)))
		self.reactorList.itemDoubleClicked.connect(self.HandleInputReactorClick)
		reacInLayout.addWidget(self.reactorList)
		self.reactorInGroup.setLayout(reacInLayout)

		self.targInGroup = QGroupBox("Input Targets", self)
		targInLayout = QVBoxLayout()
		self.targetList = QListWidget(self.targInGroup)
		for targ in self.targets:
			item = QListWidgetItem()
			self.targetList.addItem(item)
			self.targetList.setItemWidget(item, QLabel(str(targ)))
		self.targetList.itemDoubleClicked.connect(self.HandleInputTargetClick)
		targInLayout.addWidget(self.targetList)
		self.targInGroup.setLayout(targInLayout)

		self.selectedGroup = QGroupBox("Selected", self)
		selectedLayout = QVBoxLayout()
		selectReacLabel = QLabel("Selected Reactions (In Order)")
		self.selectReacList = QListWidget(self.selectedGroup)
		self.selectReacList.itemDoubleClicked.connect(self.HandleSelectReactorClick)
		selectTargDecoLabel = QLabel("Selected Target:")
		self.selectTargLabel = QLabel("")
		self.selectTargClearButton = QPushButton("Clear Sel. Target", self.selectedGroup)
		self.selectTargClearButton.clicked.connect(self.HandleClearSelTarget)
		selectedLayout.addWidget(selectReacLabel)
		selectedLayout.addWidget(self.selectReacList)
		selectedLayout.addWidget(selectTargDecoLabel)
		selectedLayout.addWidget(self.selectTargLabel)
		selectedLayout.addWidget(self.selectTargClearButton)
		self.selectedGroup.setLayout(selectedLayout)

		listLayout.addWidget(self.reactorInGroup)
		listLayout.addWidget(self.targInGroup)
		listLayout.addWidget(self.selectedGroup)
		self.layout.addWidget(self.listGroup)

	def HandleInputReactorClick(self, item):
		widget = self.reactorList.itemWidget(item)
		new_item = QListWidgetItem()
		row = self.reactorList.row(item)
		self.selectReacList.addItem(new_item)
		self.selectReacList.setItemWidget(new_item, QLabel(widget.text()))
		self.selectReactors.append(self.reactors[row])

	def HandleSelectReactorClick(self, item):
		row = self.selectReacList.row(item)
		self.selectReacList.takeItem(row)
		self.selectReactors.pop(row)
		del item
		return

	def HandleInputTargetClick(self, item):
		widget = self.targetList.itemWidget(item)
		self.selectTargLabel.setText(widget.text())
		self.selectedTarget = self.targets[self.targetList.row(item)]

	def HandleClearSelTarget(self):
		self.selectTargLabel.setText("")
		self.selectedTarget = None

	def SetInitialValues(self, chain):
		for reactor in chain.reactors:
			item = QListWidgetItem()
			self.selectReacList.addItem(item)
			self.selectReacList.setItemWidget(item, QLabel(str(reactor)))
		if chain.target is not None:
			self.selectTargLabel.setText(str(chain.target))

	def SendChain(self):
		chain = ReactorChain(self.selectReactors, self.selectedTarget)
		if chain.VerifyChain():
			self.new_chain.emit(chain)
		else:
			print("Chain failed the VerifyChain method!")
		return

class ReactorDialog(QDialog):
	new_reactor = Signal(Reactor)
	def __init__(self, parent=None, reactor=None):
		super().__init__(parent)
		self.setWindowTitle("Add A Reaction")

		QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		self.buttonBox = QDialogButtonBox(QBtn)
		self.buttonBox.accepted.connect(self.accept)
		self.buttonBox.accepted.connect(self.SendReactor)
		self.buttonBox.rejected.connect(self.reject)
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

		self.type = ReactorType.Reaction

		self.CreateReactorInputs()
		if reactor is not None:
			self.SetInitialValues(reactor)
		self.layout.addWidget(self.buttonBox)

	def CreateReactorInputs(self):
		self.typeButtonGroup = QGroupBox("Reaction/Decay Switch",self)
		buttonLayout = QHBoxLayout()
		self.rxnButton = QRadioButton("Reaction", self.typeButtonGroup)
		self.rxnButton.toggle()
		self.rxnButton.toggled.connect(self.HandleRxnSwitch)
		self.decayButton = QRadioButton("Decay", self.typeButtonGroup)
		self.decayButton.toggled.connect(self.HandleDecaySwitch)
		buttonLayout.addWidget(self.rxnButton)
		buttonLayout.addWidget(self.decayButton)
		self.typeButtonGroup.setLayout(buttonLayout)
		self.layout.addWidget(self.typeButtonGroup)

		self.inputGroup = QGroupBox("Inputs", self)
		inputLayout = QHBoxLayout()

		self.nucleiGroup = QGroupBox("Nuclei", self.inputGroup)
		nucleiLayout = QFormLayout()
		self.ztInput = QSpinBox(self.nucleiGroup)
		self.ztInput.setRange(1, 110)
		self.atInput = QSpinBox(self.nucleiGroup)
		self.atInput.setRange(1,270)
		self.zpInput = QSpinBox(self.nucleiGroup)
		self.zpInput.setRange(1, 110)
		self.apInput = QSpinBox(self.nucleiGroup)
		self.apInput.setRange(1,270)
		self.zeInput = QSpinBox(self.nucleiGroup)
		self.zeInput.setRange(1, 110)
		self.aeInput = QSpinBox(self.nucleiGroup)
		self.aeInput.setRange(1,270)
		nucleiLayout.addRow("ZT",self.ztInput)
		nucleiLayout.addRow("AT",self.atInput)
		nucleiLayout.addRow("ZP",self.zpInput)
		nucleiLayout.addRow("AP",self.apInput)
		nucleiLayout.addRow("ZE",self.zeInput)
		nucleiLayout.addRow("AE",self.aeInput)
		if self.type == ReactorType.Decay:
			self.zpInput.setEnabled(False)
			self.apInput.setEnabled(False)

			self.beamMeanInput.setEnabled(False)
			self.beamSigmaInput.setEnabled(False)
		self.nucleiGroup.setLayout(nucleiLayout)
		inputLayout.addWidget(self.nucleiGroup)

		self.beamGroup = QGroupBox("Reaction Params.", self.inputGroup)
		beamLayout = QFormLayout()
		self.beamMeanInput = QDoubleSpinBox(self.beamGroup)
		self.beamMeanInput.setDecimals(4)
		self.beamSigmaInput = QDoubleSpinBox(self.beamGroup)
		self.beamSigmaInput.setDecimals(4)
		self.exMeanInput = QDoubleSpinBox(self.beamGroup)
		self.exMeanInput.setDecimals(4)
		self.exSigmaInput = QDoubleSpinBox(self.beamGroup)
		self.exSigmaInput.setDecimals(4)
		beamLayout.addRow("Mean BeamE (MeV)", self.beamMeanInput)
		beamLayout.addRow("Sigma BeamE (MeV)", self.beamSigmaInput)
		beamLayout.addRow("Mean Ex (MeV)", self.exMeanInput)
		beamLayout.addRow("Sigma Ex (MeV)", self.exSigmaInput)
		self.beamGroup.setLayout(beamLayout)
		inputLayout.addWidget(self.beamGroup)
	
		self.inputGroup.setLayout(inputLayout)
		self.layout.addWidget(self.inputGroup)

	def HandleRxnSwitch(self):
		if self.rxnButton.isChecked() and (not self.type == ReactorType.Reaction):
			self.type = ReactorType.Reaction
			self.zpInput.setEnabled(True)
			self.apInput.setEnabled(True)
			self.beamMeanInput.setEnabled(True)
			self.beamSigmaInput.setEnabled(True)

	def HandleDecaySwitch(self):
		if self.decayButton.isChecked() and (not self.type == ReactorType.Decay):
			self.type = ReactorType.Decay
			self.zpInput.setEnabled(False)
			self.apInput.setEnabled(False)
			self.beamMeanInput.setEnabled(False)
			self.beamSigmaInput.setEnabled(False)

	def SetInitialValues(self, reactor):
		if reactor.type == ReactorType.Reaction:
			self.type = ReactorType.Reaction
			self.ztInput.setValue(reactor.nuclei[0].Z)
			self.atInput.setValue(reactor.nuclei[0].A)
			self.zpInput.setValue(reactor.nuclei[1].Z)
			self.apInput.setValue(reactor.nuclei[1].A)
			self.zeInput.setValue(reactor.nuclei[2].Z)
			self.aeInput.setValue(reactor.nuclei[2].A)
			self.beamMeanInput.setValue(reactor.beamEnergyMean)
			self.beamSigmaInput.setValue(reactor.beamEnergySigma)
			self.exMeanInput.setValue(reactor.residualExMean)
			self.exSigmaInput.setValue(reactor.residualExSigma)
			self.zpInput.setEnabled(True)
			self.apInput.setEnabled(True)
			self.beamMeanInput.setEnabled(True)
			self.beamSigmaInput.setEnabled(True)
		elif reactor.type == ReactorType.Decay:
			self.type = ReactorType.Decay
			self.ztInput.setValue(reactor.nuclei[0].Z)
			self.atInput.setValue(reactor.nuclei[0].A)
			self.zeInput.setValue(reactor.nuclei[1].Z)
			self.aeInput.setValue(reactor.nuclei[1].A)
			self.exMeanInput.setValue(reactor.residualExMean)
			self.exSigmaInput.setValue(reactor.residualExSigma)
			self.zpInput.setEnabled(False)
			self.apInput.setEnabled(False)
			self.beamMeanInput.setEnabled(False)
			self.beamSigmaInput.setEnabled(False)

	def SendReactor(self):
		reac = None
		nuclei = []
		if self.type == ReactorType.Reaction:
			nuclei.append(Nucleus(self.ztInput.value(), self.atInput.value()))
			nuclei.append(Nucleus(self.zpInput.value(), self.apInput.value()))
			nuclei.append(Nucleus(self.zeInput.value(), self.aeInput.value()))
			reac = Reactor(nuclei, self.exMeanInput.value(), self.exSigmaInput.value(), self.beamMeanInput.value(), self.beamSigmaInput.value())
		elif self.type == ReactorType.Decay:
			nuclei.append(Nucleus(self.ztInput.value(), self.atInput.value()))
			nuclei.append(Nucleus(self.zeInput.value(), self.aeInput.value()))
			reac = Reactor(nuclei, self.exMeanInput.value(), self.exSigmaInput.value())

		if reac is not None and reac.CheckNuclei():
			self.new_reactor.emit(reac)
		else:
			print("Invalid reaction!")

class TargetDialog(QDialog):
	new_target = Signal(Target)
	def __init__(self, parent=None, target=None):
		super().__init__(parent)
		self.setWindowTitle("Add A Target")

		QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		self.buttonBox = QDialogButtonBox(QBtn)
		self.buttonBox.accepted.connect(self.accept)
		self.buttonBox.accepted.connect(self.SendTarget)
		self.buttonBox.rejected.connect(self.reject)
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

		self.CreateTargetInputs()
		if target is not None:
			self.SetInitialValues(target)
		self.layout.addWidget(self.buttonBox)

	def CreateTargetInputs(self):

		self.ZInputs = []
		self.SInputs = []

		self.targGroupBox = QGroupBox("Parameters", self)
		targetLayout = QVBoxLayout()
		thickLabel = QLabel("Thickness(ug/cm^2)", self.targGroupBox)
		self.thickInput = QDoubleSpinBox(self.targGroupBox)
		self.thickInput.setRange(0, 999.0)
		self.thickInput.setDecimals(4)
		self.componentBox = QGroupBox("Components", self.targGroupBox)
		compLayout = QGridLayout()
		compLayout.addWidget(QLabel("Z", self.componentBox), 0, 1)
		compLayout.addWidget(QLabel("Stoich", self.componentBox), 0, 3)
		for i in range(3):
			compLayout.addWidget(QLabel("Component"+str(i), self.componentBox), i+1, 0)
			self.ZInputs.append(QSpinBox(self.componentBox))
			self.SInputs.append(QSpinBox(self.componentBox))
			compLayout.addWidget(self.ZInputs[i], i+1, 1)
			compLayout.addWidget(self.SInputs[i], i+1, 3)
		self.componentBox.setLayout(compLayout)
		targetLayout.addWidget(thickLabel)
		targetLayout.addWidget(self.thickInput)
		targetLayout.addWidget(self.componentBox)
		self.targGroupBox.setLayout(targetLayout)

		self.layout.addWidget(self.targGroupBox)

	def SetInitialValues(self, target):
		self.thickInput.setValue(target.thickness)
		for j in range(len(target.Z)):
			self.ZInputs[j].setValue(target.Z[j])
			self.SInputs[j].setValue(target.S[j])

	def SendTarget(self):
		targ = Target()
		
		Z = []
		S = []
		z = 0
		s = 0
		thick = self.thickInput.value()
		for i in range(3):
			z = self.ZInputs[i].value()
			s = self.SInputs[i].value()
			if z != 0 and s != 0:
				Z.append(z)
				S.append(s)

		if len(Z) != 0:
			targ.Z = Z
			targ.S = S
			targ.thickness = thick
			self.new_target.emit(targ)


def main():
	myapp = QApplication(sys.argv)
	window = RoleGUI()
	sys.exit(myapp.exec_())

if __name__ == '__main__':
	main()
