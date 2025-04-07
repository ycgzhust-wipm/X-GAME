command = xplor.command
from jCoupPot import JCoupPot
from noePot import NOEPot

import prePot

from xplor import select

from xplorPot import XplorPot
from rdcPotTools import *
from pdbTool import *
from atomAction import *
from selectTools import *
from simulationTools import *
from ivm import IVM

import protocol
import monteCarlo
protocol.initRandomSeed()   #set random seed - by time


protocol.initParams('/opt/xplor/xplor-nih-2.53/toppar/parallhdg_new.pro')

protocol.initStruct('ef1g_H.psf',erase=False)
protocol.initCoords('ef1g_H.pdb')
AtomSel("all").apply( SetProperty('segmentName', 'ALT0') )


protocol.initStruct('ef1g_H.psf',erase=False)
protocol.initCoords('ef1g_H.pdb')
AtomSel("all and not (segid ALT0 )").apply( SetProperty('segmentName', 'BLT0') )

#if more conformers are needed, duplicate the correspond proteins
#command('duplicate selection=( segid BLT0) segid="BLT1" end')


protocol.initNBond(repel=1.2)

command("""

    constraints

    inter = (segid ALT0 and resid 1:228)(all and not (segid ALT0 and resid 1:228))
    inter = (segid ALT0 and resid 290:437)(all and not (segid ALT0 and resid 290:437))
    inter = (segid ALT0 and resid 229:289)(all)
    weights * 1 end end

    """)

if xplor.p_processID==0:
  command("write psf output=complex.psf end")


init_t  = 3000
final_t = 25

cool_steps = 12000

from simulationTools import MultRamp, StaticRamp, InitialParams
rampedParams=[]

potList = PotList()
potList.add( XplorPot("BOND") )

potList.add( XplorPot("ANGL") )
rampedParams.append( MultRamp(0.4,1,"potList['ANGL'].setScale(VALUE)") )

potList.add( XplorPot("IMPR") )
rampedParams.append( MultRamp(0.4,1,"potList['IMPR'].setScale(VALUE)") )

potList.add( XplorPot("VDW") )
rampedParams.append( MultRamp(1.2,0.75,
                              "command('param nbonds repel VALUE end end')") )
rampedParams.append( MultRamp(.004,4,
                              "command('param nbonds rcon VALUE end end')") )


noe=PotList('noe')
potList.append(noe)
from noePotTools import create_NOEPot

pot = create_NOEPot('xlms',"xlms.tbl")
pot.setPotType("hard")
pot.setScale(2)       
pot.setAveType("sum")
noe.append(pot)

rampedParams.append( MultRamp(2,30, "noe.setScale( VALUE )") )


dyn  = IVM()

dyn.fix("""segid BLT0  """)
dyn.fix(""" segid ALT0 and resid 290:437""")
dyn.group(""" segid ALT0 and resid 1:228 """)

def structLoopAction(loopInfo):

    protocol.initMinimize(dyn, potList=potList)
    InitialParams( rampedParams )
    dyn.run()


    ini_timestep = 0.010
    potList["VDW"].setScale(0)
    protocol.initDynamics(dyn,
                          potList=potList,
                          bathTemp=init_t,
                          initVelocities=True,
                          stepsize=ini_timestep,
                          finalTime=10,
                          printInterval=100)
    dyn.run()


    timestep=ini_timestep
    potList["VDW"].setScale(1)

    protocol.initDynamics(dyn,
                          potList=potList,
                          bathTemp=init_t,
                          initVelocities=True,
                          stepsize=timestep,
                          finalTime=0.5,
                          printInterval=100)

    dyn.setResetCMInterval( 100 )

    AnnealIVM(initTemp =init_t,
              finalTemp=final_t,
              numSteps = 50,
              ivm=dyn,
              rampedParams = rampedParams).run()
    protocol.initMinimize(dyn)
    dyn.run()

    loopInfo.writeStructure(potList)

    pass

StructureLoop(numStructures=480,
              pdbTemplate='./Calc/Calc_STRUCTURE.pdb',
              structLoopAction=structLoopAction).run()
