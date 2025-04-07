#!/bin/bash

# step1: prepare the pdb and psf file for the calculation
xplorhome=/opt/xplor/xplor-nih-2.53/toppar
pdbname=e2
#pdbname=ub
cat <<eof > psf.inp
rtf @${xplorhome}/topallhdg_new.pro end

parameter @${xplorhome}/parallhdg_new.pro end


segment
   name=" "
   SETUP=TRUE
   chain
      @${xplorhome}/toph11.pep

coor @${pdbname}.pdb
end
end
end

delete select (name OT1 or name OT2) end

write psf output=${pdbname}.psf end

stop
eof
xplor -in psf.inp

cat <<eof > hbuild.inp
rtf @${xplorhome}/topallhdg_new.pro
end

parameter @${xplorhome}/parallhdg_new.pro
end


structure @${pdbname}.psf end

coor @${pdbname}.pdb

delete  select (name H*) end


hbuild select=(name H*) phistep=360 end
hbuild select=(name H*) phistep=5 end

flags exclude * include bonds angle impr end !
constraint fix (not name H*) end
mini powell nstep 1000 end

write coor output=${pdbname}_H.pdb end
write psf output=${pdbname}_H.psf end
stop
eof
xplor -in hbuild.inp

