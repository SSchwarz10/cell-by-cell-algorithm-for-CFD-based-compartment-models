/*******************************************************************
This UDF exports the data from the CFD simulation to individual text files using ANSYS 2019R1.
The file must be compiled and ANSYS FLUENT must run on a singel core, so that all cells are on one core and the headnode can loop over all cells/faces.
The cell thread and the relevant threads for the faces (internal, inlet, outlet and wall) are assumend to be one coherent thread and are not divided into separate threads.
The UDF creates the following text files

	- 01_single_cell.txt (filename1)
		- cell_ID_1 x y z V M1 M2 M3 M4 M5 |vel| dM1dx dM1dy dM1dz
	- 02_neighbor_cell.txt (filename2)
		- cell_ID_1 cell_ID_2 massflow A (inclusive inlet and outlet)
	- 03_inlet_massflow.txt (filename3)
		- cell_ID_1 -1 massflow A
	- 04_outlet_massflow.txt (filename4)
		- cell_ID_1 -2 massflow A
	- 05_wall_pipe.txt (filename5)
		- cell_ID_1 -3 massflow A

********************************************************************/

/*Import of pre-defined makros*/
#include "udf.h"
#include "metric.h"
#include "mem.h" 

/*Selection of data query*/
DEFINE_ON_DEMAND(neighbours_to_file)
	{
	
	/*****************************************************************************/  	
	/********************* Declaration of global variables ***********************/
	/*****************************************************************************/
	
		/*declaration of domain*/
		Domain *domain1=Get_Domain(1);								

		/*ID of domain*/

		int ID_Cell_Zone = 170;

		/*definition of the cell thread and the "variable" face threads*/
		Thread *thread = Lookup_Thread(domain1, ID_Cell_Zone);
 		Thread *tf;

		/*definition of flie thread*/
		FILE *fp = NULL;
		/*definitiion of cell and face identifier*/
		cell_t c, c0, c1;
		face_t f;
		
		/*definition of filenames*/
		char filename1[] = "01_single_cell.txt";
		char filename2[] = "02_neighbor_cell.txt";
		char filename3[] = "03_inlet_massflow.txt";
		char filename4[] = "04_outlet_massflow.txt";
		char filename5[] = "05_wall_pipe.txt";

		
		int i;
 		int n = 0;
		
		/*ID of face thread for inlet*/
		int ID_IN = 46;
		
		/*ID of face thread for outlet*/
		int ID_OUT = 33;										
		
		/*ID of face thread for interior faces*/
		int ID_INTERIOR = 171;

		/*ID of face thread for pipe wall faces*/
		int ID_WALL_PIPE = 175;

		/*definition of velocity components*/
		real u, v, w;

		/*definition of massflowrate over a face, the area vector and the area of the face*/
		real massflowrate;
		real A[ND_ND];
		real area;
		
		/*definition of the position vector for a cell*/
		real x[ND_ND];

		/*definition of the cell volume and the velocity magnitude*/
		real V;
		real vel;

		/*definition of the different moments and the gradient components of the first moment*/
		real M1, M2, M3, M4, M5, variance;
		real dM1dx, dM1dy, dM1dz;
		
		/*allocate the storage for the reconstructed and normal gradient for the first UDS (index 0, first moment = mean age)*/
		Alloc_Storage_Vars(domain1, SV_UDSI_RG(0), SV_UDSI_G(0), SV_NULL);
		Scalar_Reconstruction(domain1, SV_UDS_I(0), -1, SV_UDSI_RG(0), NULL);
		Scalar_Derivatives(domain1, SV_UDS_I(0), -1, SV_UDSI_G(0), SV_UDSI_RG(0), NULL);

									
	/*****************************************************************************/		
	/********************************* Main Part *********************************/	
	/*****************************************************************************/	
		
		/****************************************/
		/*01_single_cell.txt*/
		/****************************************/
		
		if ((fp = fopen(filename1, "w"))==NULL)
       			{
			Message("\n Warning: Unable to open %s for writing\n",filename1);
			}
 		else
       			{
			Message("\nWriting Data to %s...",filename1);
			/*headline within the text file*/
			fprintf(fp, "c0	x [m]	y [m]	z [m]	V [m^3]	M1 [s]	M2[s^2]	M3 [s^3]	M4 [s^4]	M5 [s^5]	vel [m/s]	dM1dx [s/m]	dM1dy [s/m]	dM1dz [s/m]\n");
					/*Loop over all cells from cell thread "thread"*/
 					begin_c_loop(c,thread)							
						{
							/*call of the volume and the position vector*/
							V = C_VOLUME(c,thread);
							C_CENTROID(x, c, thread);

							/*call of the velocity components*/
							u = C_U(c, thread);
							v = C_V(c, thread);
							w = C_W(c, thread);

							/*velocity magnitude*/
							vel = sqrt(pow(u, 2) + pow(v, 2) + pow(w, 2));

							/*call of the differents moments (UDS) and the gradient compoents for the first moment*/
							M1 = C_UDSI(c, thread, 0);
							M2 = C_UDSI(c, thread, 1);
							M3 = C_UDSI(c, thread, 2);
							M4 = C_UDSI(c, thread, 3);
							M5 = C_UDSI(c, thread, 4);
							dM1dx = C_UDSI_G(c,thread,0)[0];
							dM1dy = C_UDSI_G(c,thread,0)[1];
							dM1dz = C_UDSI_G(c,thread,0)[2];
							/*writing the data for every cell*/
							fprintf(fp, "%d 	%g	%g	%g	%g	%g	%g	%g	%g	%g	%g	%g	%g	%g\n", c, x[0], x[1], x[2], V, M1, M2, M3, M4, M5, vel, dM1dx, dM1dy, dM1dz);
						}			
 						end_c_loop(c,thread)
			/*closing of the file*/
			fclose(fp); 							
 			Message("Done\n");
			}		
		
			/***************************************************************************************************/
			/*02_neighbor_cell.txt*/
			/***************************************************************************************************/

			if ((fp = fopen(filename2, "w")) == NULL)
			{
				Message("\n Warning: Unable to open %s for writing\n", filename2);
			}
			else
			{
				Message("\nWriting Data to %s...", filename2);
				/*headline within the text file*/
				fprintf(fp, "c0	c1	massflow [kg/s]	A [m^2]\n");

				/*loop over all face threads within the domain*/
				thread_loop_f(tf, domain1)
				{
					/*loop over all faces within the face thread "tf"*/
					begin_f_loop(f, tf)
					{
						/*get the IDs of the neighbored cells*/
						c0 = F_C0(f, tf);
						c1 = F_C1(f, tf);

						if (THREAD_ID(tf) == ID_INTERIOR)
						{
							/*get the massflowrate, the area vector and the magnitude of the area vector*/
							massflowrate = F_FLUX(f, tf);
							F_AREA(A, f, tf);
							area = NV_MAG(A);
							fprintf(fp, "%d	%d	%g %g \n", c0, c1, massflowrate, area);
						}						
					}
					end_f_loop(f, tf)
						
				}
				/*closing the file*/
				fclose(fp);
				Message("Done\n");
			}

			/****************************************************************************************/
			/*03_inlet_massflow.txt*/
			/****************************************************************************************/

			if ((fp = fopen(filename3, "w")) == NULL)
			{
				Message("\n Warning: Unable to open %s for writing\n", filename3);
			}
			else
			{
				Message("\nWriting Data to %s...", filename3);
				/*headline within the text file*/
				fprintf(fp, "c0	c1	massflow [kg/s]	A [m^2]\n");


				/*loop over all face threads within the domain*/
				thread_loop_f(tf, domain1)
				{
					/*loop over all faces within the face thread "tf"*/
					begin_f_loop(f, tf)
					{
						/*get the IDs of the neighbored cells*/
						c0 = F_C0(f, tf);
						c1 = F_C1(f, tf);

						if (THREAD_ID(tf) == ID_IN)
						{
							/*get the massflowrate, the area vector and the magnitude of the area vector*/
							massflowrate = F_FLUX(f, tf);
							F_AREA(A, f, tf);
							area = NV_MAG(A);
							fprintf(fp, "%d	%d	%g %g \n", c0, -1, massflowrate, area);
						}
					}
					end_f_loop(f, tf)
				}
				/*closing the file*/
				fclose(fp);
				Message("Done\n");
			}

			/****************************************************************************************/
			/*04_outlet_massflow.txt*/
			/****************************************************************************************/

			if ((fp = fopen(filename4, "w")) == NULL)
			{
				Message("\n Warning: Unable to open %s for writing\n", filename4);
			}
			else
			{
				Message("\nWriting Data to %s...", filename4);
				/*headline within the text file*/
				fprintf(fp, "c0	c1	massflow [kg/s]	A [m^2]\n");

				/*loop over all face threads within the domain*/
				thread_loop_f(tf, domain1)
				{
					/*loop over all faces within the face thread "tf"*/
					begin_f_loop(f, tf)
					{
						/*get the IDs of the neighbored cells*/
						c0 = F_C0(f, tf);
						c1 = F_C1(f, tf);

						if (THREAD_ID(tf) == ID_OUT)
						{
							/*get the massflowrate, the area vector and the magnitude of the area vector*/
							massflowrate = F_FLUX(f, tf);
							F_AREA(A, f, tf);
							area = NV_MAG(A);

							fprintf(fp, "%d	%d	%g %g \n", c0, -2, massflowrate, area);
						}
					}
					end_f_loop(f, tf)
				}
				/*closing the file */
				fclose(fp);
				Message("Done\n");
			}

			if ((fp = fopen(filename5, "w")) == NULL)
			{
				Message("\n Warning: Unable to open %s for writing\n", filename4);
			}
			else
			{
				Message("\nWriting Data to %s...", filename5);
				/*headline within the text file*/
				fprintf(fp, "c0	c1	massflow [kg/s]	A [m^2]\n");


				/*loop over all face threads within the domain*/
				thread_loop_f(tf, domain1)
				{
					/*loop over all faces within the face thread "tf"*/
					begin_f_loop(f, tf)
					{
						/*get the IDs of the neighbored cells*/
						c0 = F_C0(f, tf);
						c1 = F_C1(f, tf);

						if (THREAD_ID(tf) == ID_WALL_PIPE)
						{
							/*get the massflowrate, the area vector and the magnitude of the area vector*/
							massflowrate = F_FLUX(f, tf);
							F_AREA(A, f, tf);
							area = NV_MAG(A);

							fprintf(fp, "%d	%d	%g %g \n", c0, -3, massflowrate, area);
						}
					}
					end_f_loop(f, tf)
				}
				/*closing the file */
				fclose(fp);
				Message("Done\n");
			}
	}
	
