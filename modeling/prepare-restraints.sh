# step2: prepare the restraint file from XL-MS data
rm xlms.tbl
while read a b
do
cat <<eof >> xlms.tbl
echo assign ( segid ALT0 and resid $a and name ca) (segid BLT0 and resi $b and name ca) 14.0 10.0 8.0
eof
done<xlms.txt

# The xlms.txt file defines the information for the cross-linking sites. The file contains two columns of numbers, where the first column represents the cross-linking sites on protein A, and the corresponding second column represents the cross-linking sites on protein B. For example:
#
# 8 48
# 72 6
# 92 6


# The three numbers, 14.0, 10.0, and 8.0, define the shape of the constraint function, with the last number typically being the one that needs modification. 
# The sum of the first and last numbers corresponds to the maximum arm length of the cross-linker. 
# In this example, the maximum arm length of the cross-linker (the distance between the CA atoms) is defined as 22 angstroms.
