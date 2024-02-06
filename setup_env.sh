conda env create --name lipid-sorting --file=environment.yml
git clone https://github.com/marrink-lab/TS2CG1.1.git
cd TS2CG1.1
./compile.sh 
cd ../
cp TS2CG1.1/PCG TS2CG-Setup-Pipeline 
conda activate lipid-sorting
pip install -e .
