// @author Merve Asiler

#include "KernelComputation.h"
#include "SceneManager.h"
#include "CGALUtils.h"

int main(int argc, char* argv[])
{

	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
	{
		string command_type = "kertrack";	// argv[1];
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-Star-shaped-meshes/350.off"; // argv[2];
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB_ItalianStarShapes/Complex_Models/star.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB_ItalianStarShapes/Thingi/203289.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB_ItalianStarShapes/Refinements/spiral/spiral1.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB-StarCandidates/Banana.obj";
		string shape_path = "D:/VS_Workspace/3D_Databases/DB_Kids/0001.isometry.1.off";
		//string shape_path = "D:/VS_Workspace/3D_Databases/DB_Horse/1.obj";

		// DRAW:
		// Example:	draw C:/Users/Merve/3D_DATABASES/DB_Kids/method_alexa/11to15at0.5.off
		//			draw C:/Users/Merve/3D_DATABASES/DB_Camel/camel-02.obj
		if (command_type == "draw")
			drawMeshToScene(shape_path);

		// COMPARE KERNEL RESULTS for CGAL & MDF-Ker-Plus & SDF-Ker-Plus
		else if (command_type == "experiment")
			doExperimentForPaper(shape_path);

		// COMPUTE KERNEL BY FLOOD FILL FROM A SINGLE POINT -> SDF-KER
		// Example: sdfker C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "sdfker")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY FLOOD FILL FROM A SINGLE POINT -> SDF-KER_PLUS (Improvement of SDF-KER)
		// Example: sdfker_plus C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "sdfker_plus")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY FLOOD FILL FROM A SINGLE POINT -> MDF-KER
		// Example: mdfker C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "mdfker")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY FLOOD FILL FROM A SINGLE POINT -> MDF-KER_PLUS (Improvement of MDF-KER)
		// Example: mdfker_plus C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "mdfker_plus")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY FLOOD FILL FROM A SINGLE POINT -> MDF-KER_TURBO
		// Example: mdfker-turbo C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "kertrack")
			ComputeKernel(shape_path, command_type);

		// COMPUTE KERNEL BY CGAL
		// Example: kernelByCGAL C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "kernel_by_cgal")
			ComputeKernel(shape_path, command_type);

		// FIND A KERNEL POINT MAXIMIZING A STATED COST FUNCTION by THIRD PARTY LIBRARY : SDLP
		// Example: sdlp C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Banana.obj
		else if (command_type == "findkernelpoint_SDLP")
			FindKernelPoint_SDLP(shape_path);

		// APPLY SPHERICAL PARAMETRIZE by SENDING RAYS FROM A KERNEL POINT INSIDE THE MESH
		// Example: sp C:/Users/Merve/3D_DATABASES/DB-Star-shaped-meshes/Rock_6.obj
		else if (command_type == "sp")
			SphericalParametrize(shape_path);

		// EXTRACT CONVEX HULL OF MESH
		// Example:	convexhull C:/Users/Merve/3D_DATABASES/DB_FaustRegistrations/tr_reg_000.off
		else if (command_type == "convexhull")
			computeConvexHull(shape_path);

		// UNDEFINED:
		else
			cout << "Undefined Operation!" << endl;
	}
	_CrtDumpMemoryLeaks();

	return 0;

}
