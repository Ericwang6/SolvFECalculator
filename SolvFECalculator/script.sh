if [ $(which gmx_mpi) ]; then
  gmx=gmx_mpi
else
  gmx=gmx
fi

cd em
$gmx grompp -f em.mdp -c solv.gro -p processed.top -o em.tpr -maxwarn 10
$gmx mdrun -deffnm em -v
cd ../

cd nvt
$gmx grompp -f nvt.mdp -c ../em/em.gro -p processed.top -r ../em/em.gro -o nvt.tpr -maxwarn 10
$gmx mdrun -deffnm nvt
if [ $? -ne 0 ]; then
  echo "NVT Failed";
  exit 1
fi
cd ../

cd npt
$gmx grompp -f npt.mdp -c ../nvt/nvt.gro -t ../nvt/nvt.cpt -p processed.top -r ../nvt/nvt.gro -o npt.tpr -maxwarn 10
$gmx mdrun -deffnm npt
if [ $? -ne 0 ]; then
  echo "NPT Failed";
  exit 1
fi
cd ../

cd prod
$gmx grompp -f prod.mdp -c ../npt/npt.gro -t ../npt/npt.cpt -p processed.top -r ../npt/npt.gro -o prod.tpr -maxwarn 10
$gmx mdrun -deffnm prod
if [ $? -ne 0 ]; then
  echo "Production Failed";
  exit 1
fi