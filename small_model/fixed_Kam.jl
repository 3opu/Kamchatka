#Scenario 1, modified after Koulakov
#KLU  connected to the source through a pipe [cylinder]
#BEZ having 2 chambers at surface level [0-2km] and the other one at 8-10km depth. Connected directly to the source through a pipe

using GeophysicalModelGenerator
using GMT


####################################################################
# Creating the topography and a projection
data_Topo = ImportTopo(lon = [160, 161], lat=[55.5, 56.5], file="@earth_relief_03s.grd")
#Write_Paraview(data_Topo, "Topography_Kamchatka")

proj = ProjectionPoint(Lon=160.642, Lat=56.056) #Klyuchevskoy location lon/lat
Topo_cart = Convert2CartData(data_Topo, proj)
#Write_Paraview(Topo_LaMEM, "Topo2_LaMEM")

######################################################################
# Create LaMEM setup
######################################################################

Grid = ReadLaMEM_InputFile("Kamchatka.dat")         #read the LaMEM input file of the simulation.

#nProcX, nProcY, nProcZ, xc, yc, zc, nNodeX, nNodeY, nNodeZ = GetProcessorPartitioning("ProcessorPartitioning_8cpu_2.2.2.bin")  #read the partitioning file

Topo_LaMEM = CartData(XYZGrid(Grid.xn_vec,Grid.yn_vec,0.0))
Topo_LaMEM = ProjectCartData(Topo_LaMEM, data_Topo, proj)

#In a next step we need to give each of these points a Phase number 
#(which is an integer, that indicates the type of the rock that point has), as well as the temperature (in Celcius).
Phases = ones(Int32, size(Grid.X));
Temp = zeros(Float64, size(Grid.X));

Geotherm = 20;
Temp = -Grid.Z.*Geotherm;

Temp[Temp.<20] .= 20;                               #Set temperatures not lower than 20 degrees
Temp[Temp.>1350] .=1350;                            #Mantle adiabat



#Add box for air and every layer after Dobretsov
Crust = AddBox!(Phases, Temp, Grid; xlim=(-20, 20), ylim=(-20, 20), zlim=(-30, 0), Origin=nothing, phase = ConstantPhase(1),T=LinearTemp(0, 700));
# Mantle
Mantle = AddBox!(Phases, Temp, Grid; xlim=(-20, 20), ylim=(-20, 20), zlim=(-50, -30), Origin=nothing, phase = ConstantPhase(2),T=LinearTemp(700, 1000));

######### underneath KLU
AddEllipsoid!(Phases,Temp,Grid, cen=(0, 0, -40), axes=(12, 14, 5), StrikeAngle=0, DipAngle=0, phase=ConstantPhase(4), T=ConstantTemp(1200));  #Magma source?
AddCylinder!(Phases, Temp, Grid, base=(0, 0, -40), cap=(0, 0, 10), radius=0.75, phase=ConstantPhase(4), T=ConstantTemp(1200));              #Pipe underneath KLU

######## underneath BEZ
AddCylinder!(Phases, Temp, Grid, base=(-3.5, -9.3, -40), cap=(-3.5, -9.3, -13), radius=0.75, phase=ConstantPhase(4), T=ConstantTemp(1000));             #Pipe underneath BEZ
AddEllipsoid!(Phases,Temp,Grid, cen=(-3.5, -9.3, -13), axes=(4, 6, 2), StrikeAngle=0, DipAngle=0, phase=ConstantPhase(3), T=ConstantTemp(1000));   #Chamber at 10km depth BEZ 
AddCylinder!(Phases, Temp, Grid, base=(-3.5, -9.3, -13), cap=(-3.5, -9.3, 10), radius=0.75, phase=ConstantPhase(3), T=ConstantTemp(1000)); 

ind = AboveSurface(Grid, Topo_LaMEM);
Phases[ind] .= 0;

Model3D = CartData(Grid, (Phases=Phases, Temp=Temp))
#Write_Paraview(Model3D, "model_asd")

Save_LaMEMTopography(Topo_LaMEM, "Topography.txt")
Save_LaMEMMarkersParallel(Model3D; PartitioningFile="ProcessorPartitioning_32cpu_4.2.4.bin", directory="./markers", verbose=true)
