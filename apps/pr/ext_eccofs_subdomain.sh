#!/bin/bash
# Spatially subset ROMS ECCOFS files for sub-domain applications
#
# USAGE:
#     ext_eccofs_subdomain.sh <filetype> <decimate> <number>
#         filetype = grdfile || tidfile || inifile || clmfile
#         decimate = 1 or 0 to decimate the grid, or not
#         number (optional) = part of the 3-D file name if calling
#         this script from ext_driver_multifile.sh
#
# The user selects the starting and ending xi,eta index for the RHO POINTS
# grid. Indicies for _u,_v points dimensions will be computed so
# that the extracted file has the correct dimensions for a complete ROMS file
#
# John Wilkin - 26 April 2024 / rewritten 20 May 2025
#
# To have this subset grid meet the rules for use in conjunction with a factor
# of 2 decimated grid in mixed resolution DA, then the dimensions must be:
# Lmc (coarse) = (Lmf-1)/2 => Lmf = 2*Lmc+1 => Lmf and Lpf (xi_rho) are odd
# Mmc          = (Mmf-1)/2
# but this is not required if the subset is being used solely to run a forward
# nonlinear model subdomain of ECCOFS or if the intention is to subsequently
# refine the grid or create nests.
#
# User inputs ---------------------------------------------------------- 
#
# Set the spatial subset xi,eta limits of full resolution (3-km) ECCOFS grid
# For mixed-res DA grids, both must be even numbers

# Puerto Rico subdomain covering the Bahamas for hurricane Dorian 
xmin=597
xmax=998
emin=89
emax=618

# Subdomain label becomes part of the extracted filenames
subdomain=pr3km

# End user inputs ------------------------------------------------------ 

# Process script inputs

# $1 filetype
# Set the file type to be extracted: GRID, TIDE, ROMSFILE (his,avg,clm) ...            
# filetype=gridfile - extract a subset grid file
# filetype=tidefile - extract a subset tide file
# filetype=initfile - 3-D state just one time record
# filetype=climfile - 3-D state process all time records in input
filetype=$1
case $filetype in
    grdfile|tidfile|inifile|clmfile)
        ;;
    *)
	echo "filetype must be grdfile, tidfile, inifile or clmfile"
	exit
	;;
esac

# $2 decimation flag 
# 1 (true) if the output is to be decimated
decimategrid=$2

if [ "$decimategrid" = true ] ; then
    # Find the dimensions of the decimated grid
    let "xmin = $xmin/2"
    let "xmax = $xmax/2"
    let "emin = $emin/2"
    let "emax = $emax/2"
fi

# For all grids, end index for u,v,psi differs from rho
let "xmaxu = $xmax-1" 
let "emaxv = $emax-1"

# Set the input and output file names and dimension constraints --------------
# Adapt these file inputs and output names for the application

case $filetype in
    grdfile)
	echo "Subsetting grid file"
	infile=../grid/grid_eccofs_3km_08_b7_w3.nc
	outfile=grd_$subdomain.nc
	r_constraint="-d xi_rho,$xmin,$xmax -d eta_rho,$emin,$emax"
	u_constraint="-d xi_u,$xmin,$xmaxu -d eta_u,$emin,$emax"
	v_constraint="-d xi_v,$xmin,$xmax -d eta_v,$emin,$emaxv"
	p_constraint="-d xi_psi,$xmin,$xmaxu -d eta_psi,$emin,$emaxv"
	t_constraint=" "
	var_constraint=" "
	type="GRID file"
	;;
    tidfile)
	echo "Subsetting tide file"
	infile=tides_tpxo_eccofs_3km_09_b7_w2.nc
	outfile=tid_$subdomain.nc
	# tides only have rho dimension
	r_constraint="-d xi_rho,$xmin,$xmax -d eta_rho,$emin,$emax"
	u_constraint=" "
	v_constraint=" "
	p_constraint=" "
	t_constraint=" "
	var_constraint=" "
	type="TIDES file"
	;;
    inifile)
	echo "Subsetting 3-D file for initial conditions (one time record)"
	indir=/Users/wilkin/Dropbox/_roms-db/eccofs/fwd/Run07
        # fnum from driver script chooses the eccofs initial conditions
        fnum=$3
        infile=$indir/eccofs_avg_$fnum.nc
	outfile=ini_xeccofs_avg_$fnum"_"$subdomain.nc
        echo $outfile
	r_constraint="-d xi_rho,$xmin,$xmax -d eta_rho,$emin,$emax"
	u_constraint="-d xi_u,$xmin,$xmaxu -d eta_u,$emin,$emax"
	v_constraint="-d xi_v,$xmin,$xmax -d eta_v,$emin,$emaxv"
	p_constraint="-d xi_psi,$xmin,$xmaxu -d eta_psi,$emin,$emaxv"
	tmin=0
	tmax=0
	time_var=ocean_time    
	t_constraint="-d $time_var,$tmin,$tmax"
	var_constraint=" -v temp,salt,u,v,w,ubar,vbar,zeta,AKs,AKt,spherical,f,angle,pm,pn,mask_rho,mask_u,mask_v,mask_psi "
	type="INITIAL file"
	;;
    clmfile)
	echo "Subsetting 3-D file for clm/nud conditions (all time records)"
	indir=/Users/wilkin/Dropbox/_roms-db/eccofs/fwd/Run07
        # fnum from driver script loops over eccofs outputs to build subdomain clm files
        fnum=$3
        infile=$indir/eccofs_avg_$fnum.nc
	outfile=clm_avg_$fnum"_"$subdomain.nc
	r_constraint="-d xi_rho,$xmin,$xmax -d eta_rho,$emin,$emax"
	u_constraint="-d xi_u,$xmin,$xmaxu -d eta_u,$emin,$emax"
	v_constraint="-d xi_v,$xmin,$xmax -d eta_v,$emin,$emaxv"
	# some input files might not have psi dimension .. adapt accordingly
	p_constraint=" "
        # no constraint on time dimension extracts all records
	t_constraint=" "
	var_constraint=" "
	type="CLIMATOLOGY file"
	;;
    *)
	echo "filetype not recognized"
esac
echo "Reading from ..."
echo $infile
echo "Writing to ..."   
echo $outfile

if [ "$decimategrid" = 1 ] ; then
    echo "Decimating 3-km ECCOFS to 6-km"
else
    echo "Retaining full 3-km ECCOFS resolution"
fi

# Build the query and execute - overwrite any previous output file
# Presently, $var_constraint is not set above. In future, possibly
# introduce limited set of variables to extract, which might differ
# if e.g. hourly zeta,ubar,vbar is being extracted from QCK files
# and daily temp,salt from AVG files

ncks -4 -L 1 --overwrite $t_constraint $r_constraint $u_constraint \
     $v_constraint $p_constraint $var_constraint $infile $outfile

# Document
string="Extracted from East Coast Community Ocean Forecast System (ECCOFS) 3-km forecast grid "$infile" using script ext_eccofs_subdomain.sh (John Wilkin)"
echo $string
ncatted -h -a title,global,m,c,"$string" $outfile
string="Subdomain limits in ECCOFS xi="$xmin,$xmax", eta="$emin,$emax" "
if [ "$decimategrid" = 1 ] ; then
    string+=" with factor of 2 decimation"
fi
echo $string
ncatted -h -a subdomain,global,c,c,"$string" $outfile
ncatted -h -a type,global,m,c,"$type" $outfile


