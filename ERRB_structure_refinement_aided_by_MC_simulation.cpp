#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <cstdlib>
//#include <TMinuit.h>
//#include <TROOT.h>
//#include <TMath.h>
#include <time.h>
#include <ctime>


#define randomize() srand((unsigned)time(NULL))
#define random(n) (rand()%(n))

using namespace std;

const char deoxy_PDB[250] = "XXX";
const char carboxy_int_file[230] = "XXX";
const char reducedmass[250] = "XXX";
const char formfactor[250]  = "XXX";
const char exp_curve[250] = "XXX";

const int Nq = XXX, qini = XXX; 

const int num_of_atom_protein = XXX; 
const int num_of_water = XXX; 
const int num_of_rb = XXX;
const int num_of_loop = XXX;
const int HIS_101_NE2_A = XXX;
const int HIS_101_NE2_B = HIS_101_NE2_A + (num_of_atom_protein-num_of_water)/2;
const int HIS_101_rb = XXX;
const int Fe_A = XXX;
const int Fe_B = Fe_A+(num_of_atom_protein-num_of_water)/2;
const int Fe_rb = XXX; 
const double Fe_HIS_101 = XXX; 
const int HEM_A_NC = XXX; 
const int HEM_B_NC = HEM_A_NC+(num_of_atom_protein-num_of_water)/2;
const int subunit_atom = (num_of_atom_protein-num_of_water)/2;
const int ini_helix_aa[num_of_rb/2] ={XXX};

double fine_tune_crit = XXX;
const double chi2_crit = XXX;

double chi2_scale = XXX;


const int ini_loop_aa[num_of_loop/2]={10,26,44,51,83,104,125};

const int HIS_101_rb_B = HIS_101_rb + num_of_rb/2;
const int Fe_rb_B = Fe_rb + num_of_rb/2;

const double sigNC = 1.28, sigFE = 2.0, sigA = 1.2;
double sig;
//const double p1 0.0; // scaling factor for L-J force
double p1;
//const double p2 = 2e-7; // scaling factor for chi-force
//const double p2 = 4e-5;
double p2 = 15e-5;
double p3 = 120.0; // scaling factor for sysmetric force

const double p = 1.; // scaling factor for torque

int ini_rb[num_of_rb];
int fin_rb[num_of_rb];
int ini_loop[num_of_loop];
int fin_loop[num_of_loop];
int num_in_loop[num_of_loop];

double dt = 0.001;
double min_chi2;

const double v_annealing = 0.99;
const double w_annealing = 0.99;

// fitting parameters ****************************
const int num_of_parameter = 6;

int exp_index;

double chi2;
char chi2_corrected[255];
double chi2_prev;
double chi2_current;
double target_value_prev;
double target_value_current;
double R_prev;
double R_current;
double pen_col_prev;
double pen_col_current;
double penalty_RB_dist;
double penalty_collision;
double penalty_symmetry;
double penalty_unfolding;

bool indicator = true; // rotation step : true, RBMD step : false; 
bool indicator_ReadPDB = true;
bool indicator_PositionUpdate = false;

double para[num_of_parameter] = {
186.304,
121.118,2.65633,
1.63359,
0.0,
0.0
};

double para_result[num_of_parameter];

// para[0] : x movement to reach rot. axis
// para[1] : y movement to reach rot. axis
// para[2] : angle between rot. axis and x-axis
// para[3] : angle between rot. axis and y-axis
// para[4] : rotation angle along rot. axis
// para[5] : translation movement along rot. axis

//*************************************************

// Scailing factor ********************************************************

const double fp1 = 0.0;
const double fp2 = 300.0; //original version : 175.0
const double fp3 = 3.0e-4;
const double fp4 = 0.0;
//const double penalty_factor 1e-7;
const double penalty_factor = 0.1;
const double temp_factor = 0.01;
//double sf = 1.96525e-7; // scaling factor between exp and theory
//double sf = 8.58378e-8;
double sf = 7.006e-9;
//double sf = 2.3e-7;
//*************************************************************************

// Boltzman Annealing ****************************************************

double ini_temperature = 400.;
double temperature;
double cooling_rate = 0.9;
double final_temperature = 390.;

double d_initial = 0.0; // initial bond length
double d_final;
double *chi;
double *q; 
double *intensity_final;
double *intensity_initial;
double *diff_intensity;
double *diff_intensity_prev;
double *intensity_exp;
double *intensity_theory;
double *sigma;

double UE_lj = 0.0; // total L-J potential of protein
double UE_chi = 0.0; // total chi2 potential of protein
double KE_trans = 0.0; // total translational kinetic energy of protein
double KE_rotation = 0.0; // total rotational kinetic energy of protein


// variables for PDB data & atoms********

char atom[num_of_atom_protein][10];
int num_atom[num_of_atom_protein];
char a[num_of_atom_protein][5];
char na_am[num_of_atom_protein][5];
char b[num_of_atom_protein];
int num_residue[num_of_atom_protein];
double x[num_of_atom_protein];
double y[num_of_atom_protein];
double z[num_of_atom_protein];
double occu[num_of_atom_protein];
double c[num_of_atom_protein];
char d[num_of_atom_protein][5];

double mass[num_of_atom_protein];
double x_rp[num_of_atom_protein];
double y_rp[num_of_atom_protein];
double z_rp[num_of_atom_protein];
double f_lj_x[num_of_atom_protein], f_lj_y[num_of_atom_protein], f_lj_z[num_of_atom_protein];
double f_chi_x[num_of_atom_protein], f_chi_y[num_of_atom_protein], f_chi_z[num_of_atom_protein];
double f_sym_x[num_of_atom_protein], f_sym_y[num_of_atom_protein], f_sym_z[num_of_atom_protein];
double f_tot_x[num_of_atom_protein], f_tot_y[num_of_atom_protein], f_tot_z[num_of_atom_protein];
double form[num_of_atom_protein][Nq-qini+1];
double uwcrx[num_of_atom_protein], uwcry[num_of_atom_protein], uwcrz[num_of_atom_protein];
//*****************************************



// variables for RB ************************
double com_x[num_of_rb];
double com_y[num_of_rb];
double com_z[num_of_rb];
double mass_rb[num_of_rb];
double v_x[num_of_rb];
double v_y[num_of_rb];
double v_z[num_of_rb];
double w_x[num_of_rb];
double w_y[num_of_rb];
double w_z[num_of_rb];
double f_rb_tot_x[num_of_rb], f_rb_tot_y[num_of_rb], f_rb_tot_z[num_of_rb];
double torque_x[num_of_rb], torque_y[num_of_rb], torque_z[num_of_rb];
double moment_of_inertia[num_of_rb][3][3];
double moment_inverse[num_of_rb][3][3];
double Iwx[num_of_rb], Iwy[num_of_rb], Iwz[num_of_rb];
double Ax[num_of_rb], Ay[num_of_rb], Az[num_of_rb];
double unit_w_x[num_of_rb], unit_w_y[num_of_rb], unit_w_z[num_of_rb];
//*******************************************

// variables for Rotation ****************
double alpha, beta, cos_gamma, theta, phi;
double delta;
double trans_factor;

// functions ********************************
void ReadPDB();
void MakePDBMovie();
void MakePDB();
void MoveToRotAxis();
void CoordiTransformToRotAxis();
void Rotation();
void CoordiTransfromBackToOrigin();
void TranslationAlongRotAxis();
void GoBackToOrigin();
void CrysolCal();
void ReadIniPDBFile();
void CalChi2();
//void mini_chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *para, Int_t iflg);
void MakeFile();
//index for loop and rigid body
void MakeIndex();

// additional functions for RBMD

void InitialBondCheck();
void ReadFormfactor();
void ReadMass();
void CalCom();
void CalRelativePosition();
void CalForce();
void Translation();
void Rotation_Mb();
void PositionUpdate();
void AArearrange();
void Looprearrange();

double BondCheck();

double TargetFunction();
//*******************************************

bool BondBreak();
bool BondBreakPeptide();
void BringTerminus();
void ConstrainMove();
double target_dist = 3.8;
bool termi_break;
bool inter_break;
int num_break_rb[num_of_rb];
int num_of_break;
void SelectReplica();
double CalRfactor();

char replica_line[500];

double least_chi2_table[1000];
int least_chi2_index;

double minimum;
bool first_time;
int main(){
	ofstream out;
	int i;
	int success_step = 0;
	int repeat_step = 0;
//	system("mv 4SDH_backup.pdb ./4SDH_17H2O.pdb");
	system("rm ./result/*");
	system("rm ./result_pdb/*");
	system("rm movie.pdb");
	system("rm ./new*");
	system("rm ./candidate/*");
	system("rm -r replica");
	system("cp -r replica_pool ./replica");
	system("mkdir minimum_chi2_table");
	system("mkdir ./result_pdb");

	system("chmod 777 crysol26.l86");
	system("chmod 777 polychromatic_correction_harry");
//	system("cp ./4SDH_17H2O.pdb ./4SDH_backup.pdb");
	
	MakeIndex();
	
	intensity_final = new double[Nq+1];
	intensity_initial = new double[Nq+1];
	diff_intensity = new double[Nq+1];
	diff_intensity_prev = new double[Nq+1];
	q = new double[Nq+1];
	intensity_exp = new double[Nq-1];
	sigma = new double[Nq-1];
	intensity_theory = new double[Nq-1];
	double *intensity_theory_prev, *double_difference;
			intensity_theory_prev = new double[Nq-1];
            double_difference = new double[Nq+1];
	chi = new double[Nq-qini+1];

	for(i=0 ; i<=Nq-qini ; i++)
		chi[i] = 0.0;         


//	indicator = false;
//	indicator_ReadPDB = false;

	ReadPDB();
	MakePDB();
	InitialBondCheck();
	ReadIniPDBFile();
	ReadFormfactor();
	ReadMass();
	CalCom();
	CalRelativePosition();
	temperature = ini_temperature;


	randomize();
	cout<<"randomize()"<<endl;

	CrysolCal();
	target_value_prev = TargetFunction();
	chi2_prev = chi2;
	
	system("mv new.pdb ./prev.pdb");
	system("rm new*");


	static int count_step = 0;

	
	bool bond_break1;
	bool bond_break2;
	bool bond_break;

	indicator_ReadPDB = false;

// Initial Selection of Replica
	SelectReplica();
	CrysolCal();
	target_value_prev = TargetFunction();
	pen_col_prev = penalty_collision;
	chi2_prev = chi2;
	bond_break1 = BondBreak();
	bond_break2 = BondBreakPeptide();
	bond_break = bond_break1 + bond_break2;

	system("rm new0*");

//****************************************
	for(int j=0 ; j<10 ; j++){
		indicator = false;
		indicator_ReadPDB = false;

//		temperature = ini_temperature-(50*count_step);
		count_step++;

		bool update = true;
		double pre_p2;
		bool shock = false;


		update = true;
		shock = false;
		
//  Rolling the bead **********************************************************
		int num_of_try = 0;
		int trap_limit = 30;

//*****************************************************************************

		bool fine_tuning_mode;
		int fine_tuning_num;
		bool log_book = true;			

		while(temperature >= final_temperature){
			
				dt = 0.001;
				if(repeat_step==1)fine_tuning_mode = false;
				//if(repeat_step==1)fine_tuning_mode = true;

				static int trap_num=0;

					p2 = (random(60000)/10000.)*1e-6;
//					p2 = ((60000)/10000.)*1e-6;


					if(min_chi2<10.0) p1=5e6;
					else if(min_chi2<100.0) p1=5e5;
					else if(min_chi2<1000.0) p1=5e4;
					else p1=5e3;

					if(trap_num>trap_limit){
						num_of_try++;
						p2=6.0e-6*num_of_try;

/*						if(num_of_try!=1){
							out.open("./result/target_value_real.dat",ios::app);
							out<<"num_of_try : "<<num_of_try<<'\t'<<"trap.pdb is loaded."<<endl;
							out.close();
							system("cp ./trap_lib/trap.pdb ./new.pdb");
							ReadPDB();
							CalCom();
							CalRelativePosition();
						}
*/
						if(first_time==false){
								SelectReplica();
								fine_tuning_mode = false;
								//fine_tuning_mode = true;
								CrysolCal();
								target_value_current = TargetFunction();
								pen_col_current = penalty_collision;
								chi2_current = chi2;
								num_of_try = 0;
								repeat_step = success_step = 1;
								trap_num = 0;

								bond_break1 = BondBreak();
								bond_break2 = BondBreakPeptide();
								bond_break = bond_break1 + bond_break2;

								goto trap_case;
						}
						first_time = false;
						min_chi2 = 1e10;
					}

					MakePDB();
					
					for(int kk=0;kk<num_of_atom_protein; kk++){
						if(y[kk]<=-100 || y[kk]>=1000){
							out.open("result/collision.dat", ios::app);
							out<<"y : "<<y[kk]<<endl;
							cout<<"y : "<<y[kk]<<endl;
							out.close();
							system("rm new.pdb");
							system("cp prev.pdb ./new.pdb");
							goto strange_pdb;
						}
						if(z[kk]<=-100 || z[kk]>=1000){
							out.open("result/collision.dat", ios::app);
							out<<"z : "<<z[kk]<<endl;
							cout<<"z : "<<z[kk]<<endl;
							system("rm new.pdb");
							system("cp prev.pdb ./new.pdb");
							out.close();
							goto strange_pdb;
						}
					}
					
					CrysolCal();
					CalChi2();
					CalForce();
			
					system("rm new*");

					Translation();
					Rotation_Mb();			
					PositionUpdate();
			
/*					double v_rb_avr_tot;
					v_rb_avr_tot=0.0;
					out.open("v_rb.dat");           
					out<<p2<<endl;
					for(int i=0 ; i<num_of_rb ; i++){
                        out<<v_x[i]<<'\t'<<v_y[i]<<'\t'<<v_z[i]<<endl;
						v_rb_avr_tot+=sqrt(v_x[i]*v_x[i]+v_y[i]*v_y[i]+v_z[i]*v_z[i]);
						
					}
					out<<v_rb_avr_tot/num_of_rb<<endl;
					out.close();exit(1);
					*/
					for(int kk=0;kk<num_of_atom_protein; kk++){
						if(y[kk]<=-100 || y[kk]>=1000){
							out.open("result/collision.dat", ios::app);
							out<<"y : "<<y[kk]<<endl;
							cout<<"y : "<<y[kk]<<endl;
							out.close();
							system("rm new.pdb");
							system("cp prev.pdb ./new.pdb");
							goto strange_pdb;
						}
						if(z[kk]<=-100 || z[kk]>=1000){
							out.open("result/collision.dat", ios::app);
							out<<"z : "<<z[kk]<<endl;
							cout<<"z : "<<z[kk]<<endl;
							system("rm new.pdb");
							system("cp prev.pdb ./new.pdb");
							out.close();
							goto strange_pdb;
						}
					}

					
					MakePDB();
					ReadPDB();

					bond_break1 = BondBreak();
					bond_break2 = BondBreakPeptide();
					bond_break = bond_break1 + bond_break2;


					CrysolCal();
					target_value_current = TargetFunction();
					chi2_current = chi2;

					
					if(min_chi2<fine_tune_crit)
                                fine_tuning_mode = true;
					else
                                fine_tuning_mode = false;
	
					/************** Fine tuning mode **************************************************************************/
					
					if(fine_tuning_mode&&(trap_num!=0)){
						cout<<"fine_tuning_mode"<<endl;	
	
						if((target_value_current >= target_value_prev) || (bond_break==true)){
							system("rm new0*");
							system("cp prev.pdb ./new.pdb");
							ReadPDB();
							CalCom();
							CalRelativePosition();

							for(int i=0 ; i<=Nq ; i++)
								double_difference[i] = diff_intensity[i] - diff_intensity_prev[i];

							double min, dt_min;
							double temp[121];

							for(int i=1 ; i<=1000 ; i++){
									dt = i/100000.0;
						
									for(int k=0 ; k<=Nq ; k++){
										temp[k] = double_difference[k];
										temp[k] = temp[k]*dt/0.001*0.904;
										diff_intensity[k] = diff_intensity_prev[k] + temp[k];

									}
									chi2=0.0;
									for(int k=qini ; k<=100 ; k++){
											chi2 += pow((intensity_exp[exp_index+(k-qini)]-diff_intensity[k]),2)/pow(sigma[exp_index+(k-qini)],2);
									}
									chi2 *= chi2_scale;

									if(i==5){
										min = chi2;
										dt_min = dt;
									}
									else{
										if(min>chi2){
											min = chi2;
											dt_min = dt;
										}
									}
							}

							dt = dt_min;

							cout<<"dt_min : "<<dt_min<<endl;

							Translation();
							Rotation_Mb();
							PositionUpdate();

							for(int kk=0;kk<num_of_atom_protein; kk++){
								if(y[kk]<=-100 || y[kk]>=1000){
									out.open("result/collision.dat", ios::app);
									out<<"y : "<<y[kk]<<endl;
									cout<<"y : "<<y[kk]<<endl;
									out.close();
									system("rm new.pdb");
									system("cp prev.pdb ./new.pdb");
									goto strange_pdb;
								}
								if(z[kk]<=-100 || z[kk]>=1000){
									out.open("result/collision.dat", ios::app);
									out<<"z : "<<z[kk]<<endl;
									cout<<"z : "<<z[kk]<<endl;
									system("rm new.pdb");
									system("cp prev.pdb ./new.pdb");
									out.close();
									goto strange_pdb;
								}
							}

							MakePDB();
							ReadPDB();

							bond_break1 = BondBreak();
							bond_break2 = BondBreakPeptide();
							bond_break = bond_break1 + bond_break2;

							CrysolCal();
							target_value_current = TargetFunction();
							chi2_current = chi2;
						}
					}
					
					R_current = CalRfactor();
					pen_col_current = penalty_collision;
					if(trap_num>trap_limit){
						trap_num = 0;
						goto trap_case;
					}
					
					/************** Fine Tuning Mode End *********************************************/
					repeat_step++;
					success_step++;

					ReadPDB();

					bond_break1 = BondBreak();
					bond_break2 = BondBreakPeptide();
					bond_break = bond_break1 + bond_break2;

					static int model_num = 1;

					if((target_value_current >= target_value_prev) || (bond_break==true)){
						cout<<"update failure"<<endl;

						success_step--;

						
						if((chi2<chi2_crit)&&(penalty_collision==0.0)&&(bond_break==false)){

								MakePDBMovie();
								char temp[100];
								char cmd[100] = "cp new.pdb ";
								sprintf(temp, "./result_pdb/%04d.pdb",model_num++);
								strcat(cmd, temp);	
								system(cmd);

								out.open("./result/target_value_real.dat", ios::app);
								sprintf(chi2_corrected, "%06.4f", chi2_current);
								out<<"Model "<<model_num-1<<'\t'<<target_value_current<<'\t'<<chi2_corrected<<"\tp2 :"<<p2<<"\tR-factor : "<<R_current<<"\tpen_col : "<<pen_col_current<<"\tfine tuning : "<<fine_tuning_mode<<"\tfirst time : "<<first_time<<endl;
								out.close();
								log_book = false;
							

								char name[100];
								sprintf(name,"candidate_%06.4f_chi2.pdb",chi2);
								char cd[100] = "cp new.pdb ./candidate/";
								strcat(cd,name);
								system(cd);
								char dc[100] = "cp diff.int ./candidate/";
								char poly_name[50];
								sprintf(poly_name,"candidate_%06.4f_diff.int",chi2);
								strcat(dc,poly_name);
								system(dc);
								ofstream out;
								out.open("./candidate/log_book.log",ios::app);
								sprintf(chi2_corrected, "%06.4f", chi2);
								out<<replica_line<<" : chi2 : "<<chi2_corrected<<endl;
								out.close();
						}
						

					
strange_pdb:
						system("rm new.pdb");
						system("mv prev.pdb ./new.pdb");
						system("cp new.pdb ./trap_lib/trap.pdb");
						ReadPDB();
						CalCom();
						CalRelativePosition();

						for(i=0 ; i<num_of_rb ; i++){
      						v_x[i] = 0.0;
      						v_y[i] = 0.0;
      						v_z[i] = 0.0;
      						w_x[i] = 0.0;
      						w_y[i] = 0.0;
      						w_z[i] = 0.0;
						
						}

						indicator_PositionUpdate = true;
						PositionUpdate();

						MakePDB();
						ReadPDB();

						bond_break1 = BondBreak();
						bond_break2 = BondBreakPeptide();
						bond_break = bond_break1 + bond_break2;

				
						target_value_current = target_value_prev;
						chi2_current = chi2_prev;
						R_current = R_prev;
						
						chi2 = chi2_prev;
						penalty_collision = pen_col_prev;

						trap_num++;

						
						pen_col_current = pen_col_prev;
					}

					
					else{


						if(min_chi2<fine_tune_crit){
							fine_tuning_mode = true;
							fine_tuning_num = 0;
						}
						else 
							fine_tuning_mode = false;

						least_chi2_table[least_chi2_index++] = chi2;

						system("rm ./trap_lib/*");
						num_of_try = 0;

trap_case:

						
						for(int k=0 ; k<=Nq ; k++)
							diff_intensity_prev[k] = diff_intensity[k];

						trap_num=0;
						MakePDBMovie();
						char temp[100];
						char cmd[100] = "cp new.pdb ";
						sprintf(temp, "./result_pdb/%04d.pdb",model_num++);
						strcat(cmd, temp);	
						system(cmd);
						
						if((chi2<=chi2_crit)&&(penalty_collision==0.0)&&(bond_break==false)){
							char name[100];
							sprintf(name,"candidate_%06.4f_chi2.pdb",chi2);
							char cd[100] = "cp new.pdb ./candidate/";
							strcat(cd,name);
							system(cd);
							char dc[100] = "cp diff.int ./candidate/";
							char poly_name[100];
							sprintf(poly_name,"candidate_%06.4f_diff.int",chi2);
							strcat(dc,poly_name);
							system(dc);
							ofstream out;
							out.open("./candidate/log_book.log",ios::app);
							sprintf(chi2_corrected, "%06.4f", chi2);
							out<<replica_line<<" : chi2 : "<<chi2_corrected<<endl;
							out.close();
						}
						
					}	

					indicator_PositionUpdate = false;

					system("rm prev.pdb");
					system("cp ./new.pdb ./prev.pdb");
					target_value_prev = target_value_current;
					chi2_prev = chi2_current;
					R_prev = R_current;
					pen_col_prev = pen_col_current;
					ReadPDB();
					

					if(log_book == true)
					{
						sprintf(chi2_corrected, "%06.4f", chi2_prev);
						out.open("./result/target_value_real.dat", ios::app);
						out<<"Model "<<model_num-1<<'\t'<<target_value_prev<<'\t'<<chi2_corrected<<"\tp2 :"<<p2<<"\tR-factor : "<<R_prev<<"\tpen_col : "<<pen_col_prev<<"\tfine tuning : "<<fine_tuning_mode<<"\tfirst time : "<<first_time<<endl;
						if(bond_break)out<<"Bond Break"<<endl;
						out.close();
					}
					log_book = true;

					system("rm new*");
					
					
					if( (success_step==10000) || (repeat_step==10000) ){
						SelectReplica();
						CrysolCal();
						target_value_current = TargetFunction();
						num_of_try = 0;
						pen_col_current = penalty_collision;
						repeat_step = success_step = 1;

						bond_break1 = BondBreak();
						bond_break2 = BondBreakPeptide();
						bond_break = bond_break1 + bond_break2;

						fine_tuning_mode = false; 

						cout << "you are in unwanted place" << endl;
						goto trap_case;

					}
					
					
					if(target_value_prev<minimum)
						minimum = target_value_prev;

					if(chi2_prev<min_chi2)
						min_chi2 = chi2_prev;
						
		}       

		ofstream out;
		out.open("./result/final_temperature.dat");
		out<<(temperature*10/9)<<"K"<<endl;
		out.close();
		MakePDB();


	}


	system("mv new.pdb ./HbI_late_PC_result.pdb");		
	delete intensity_final;
	delete intensity_initial;
	delete diff_intensity;
	delete q;
	delete intensity_exp;
	delete sigma;
	delete intensity_theory;
	delete double_difference;
	delete intensity_theory_prev;

	return 0;
}

void MakeIndex(){
	cout<<"MakeIndex()"<<endl;
	indicator_ReadPDB = true;
	ReadPDB();
	
	ofstream out;
	char file_name[500];
	system("mkdir ./helix_read/");
	sprintf(file_name,"helix_read/%s_loop_helix.dat",deoxy_PDB);
	out.open(file_name);
	/****************initial rb*********************/
	out<<"ini_rb_index"<<endl;
	for(int j=0; j<num_of_rb/2; j++){
		for(int i=0 ; i<subunit_atom ; i++){
			if(num_residue[i]==ini_helix_aa[j]){
				ini_rb[j]=i;
				out<<i<<",";
				
				break;
			}			
		}
	}
	for(int j=0; j<num_of_rb/2; j++){
		for(int i=subunit_atom ; i<num_of_atom_protein ; i++){
			if(num_residue[i]==ini_helix_aa[j]){
				ini_rb[j+num_of_rb/2]=i;
				out<<i<<",";
				break;
			}			
		}
	}
	

	/****************end rb*********************/
	out<<"\n\n"<<"end_rb_index"<<endl;
	for(int j=0; j<num_of_rb-1; j++){
		fin_rb[j]=ini_rb[j+1]-1;
		out<<ini_rb[j+1]-1<<",";
	}
	fin_rb[num_of_rb-1]=num_of_atom_protein-num_of_water-1;
	out<<fin_rb[num_of_rb-1]<<endl;
			
	/****************initial loop*********************/
	out<<"\n\n"<<"ini_loop_index"<<endl;
	for(int j=0; j<num_of_loop/2; j++){
		for(int i=0 ; i<subunit_atom ; i++){
			if(num_residue[i]==ini_loop_aa[j]){
				ini_loop[j]=i;
				out<<i<<",";
				break;
			}			
		}
	}
	for(int j=0; j<num_of_loop/2; j++){
		for(int i=subunit_atom ; i<num_of_atom_protein ; i++){
			if(num_residue[i]==ini_loop_aa[j]){
				ini_loop[j+num_of_loop/2]=i;
				out<<i<<",";
				break;
			}			
		}
	}
	
	/****************end loop*********************/
	out<<"\n\n"<<"end_loop_index"<<endl;
	for(int i=0; i<num_of_loop/2; i++){
		fin_loop[i]=fin_rb[i];
		out<<fin_rb[i]<<",";
	}
	for(int i=num_of_loop/2; i<num_of_loop; i++){
		fin_loop[i]=fin_rb[i+2];
		out<<fin_rb[i+2]<<",";
	}
	out<<endl;
	out.close();
	
	/**********loop num of aa*****************/
	for(int i=0; i<num_of_loop/2; i++){
		num_in_loop[i]=ini_helix_aa[i+1]-ini_loop_aa[i];
		num_in_loop[i+num_of_loop/2]=num_in_loop[i];
	}
	
	
}


void ReadPDB(){
	cout<<"ReadPDB()"<<endl;
	ifstream in;

	if(indicator_ReadPDB==false)
		in.open("new.pdb");

	else
		in.open(deoxy_PDB);

	char *divider;
	char temp[20];
	int i=0;

	for(i=0 ; i<num_of_atom_protein ; i++){
		in>>atom[i]>>num_atom[i]>>temp;
	
		if(strlen(temp) >= 4){
			divider = temp+3;
			strncpy(a[i], temp, 3);
			a[i][3] = '\0';
			strcpy(na_am[i], divider);
		}

		else{
			strcpy(a[i], temp);
			in>>na_am[i];
		}
		
		in>>b[i]>>num_residue[i]>>x[i]>>y[i]>>z[i]>>occu[i]>>c[i]>>d[i];
	}

	in.close();
}

void MakePDB(){
	cout<<"MakePDB()"<<endl;
	FILE *stream;
	stream = fopen("new.pdb", "w+");
		for(int j=0 ; j<num_of_atom_protein ; j++){
			 if(!strcmp(a[j],"FE"))
                                fprintf(stream, "%-6s%5d %-3s %4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s  \n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
                        else if(j==num_of_atom_protein-1)
                                fprintf(stream, "%-6s%5d  %-3s%4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
                        else
                        fprintf(stream, "%-6s%5d  %-3s%4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s  \n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
		}

	fclose(stream);
}

void MakePDBMovie(){
        cout<<"MakePDB()"<<endl;
        FILE *stream;
		static int num_of_model =1;
        stream = fopen("movie.pdb", "a+");
		fprintf(stream, "MODEL%9d\n",num_of_model++);
                for(int j=0 ; j<num_of_atom_protein ; j++){
                         if(!strcmp(a[j],"FE"))
                                fprintf(stream, "%-6s%5d %-3s %4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s  \n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
                        else if(j==num_of_atom_protein-1)
                                fprintf(stream, "%-6s%5d  %-3s%4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
                        else
                        fprintf(stream, "%-6s%5d  %-3s%4s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s  \n",atom[j], num_atom[j], a[j], na_am[j], b[j], num_residue[j],x[j], y[j],z[j] ,occu[j], c[j], d[j]);
                }
		fprintf(stream,"TER\nENDMDL\n");
        fclose(stream);
}

void CrysolCal(){
	cout<<"CrysolCal()"<<endl;
	system("./crysol26.l86 <0334_003.dat");
	ifstream in;
	in.open("new00.int");
	char temp[500];
	double trash[3];
	in.getline(temp,500,'\n');
	
	for(int i=0 ; i<=Nq ; i++){
		in>>q[i]>>intensity_final[i]>>trash[0]>>trash[1]>>trash[2];
	}

	in.close();
	cout<<"CrysolCal() finished"<<endl;
}

void ReadIniPDBFile(){
	cout<<"ReadIniPDBFile()"<<endl;
	ifstream in;
	in.open(carboxy_int_file);

	char temp[500];
	double trash[3];
	in.getline(temp,500,'\n');
	
	for(int i=0 ; i<=Nq ; i++){
		in>>q[i]>>intensity_initial[i]>>trash[0]>>trash[1]>>trash[2];
	}
	in.close();
}

void CalChi2(){
	cout<<"CalChi2()"<<endl;
	int i;

	ofstream out;
	out.open("diff.int");
	for(i=0 ; i<=Nq ; i++){
		diff_intensity[i] = intensity_final[i] - intensity_initial[i];
		diff_intensity[i] *= sf;
		//intensity_theory[i] = diff_intensity[i];
		out<<q[i]<<'\t'<<diff_intensity[i]<<endl;
	}
	out.close(); // calculate scattering intensity difference between two model 
	
	double temp;

	ifstream in;
	in.open(exp_curve);
	if(in==NULL)
	{
		cerr<<"exp_curve file is not open"<<endl;
		exit(1);
	}
	for(i=0 ; i<Nq-1 ; i++){
		in>>temp>>intensity_exp[i]>>sigma[i];
		if(temp == q[qini])
			exp_index = i;			
	}
	in.close();


	if(indicator == false){
		for(i=qini ; i<=100 ; i++){
			if(min_chi2<10.0)
					chi[i-qini] = (intensity_exp[exp_index+(i-qini)]-diff_intensity[i])/pow(sigma[exp_index+(i-qini)]/sqrt(1000.0),2);
			else if(min_chi2<100.0)
					chi[i-qini] = (intensity_exp[exp_index+(i-qini)]-diff_intensity[i])/pow(sigma[exp_index+(i-qini)]/sqrt(100.0),2);
			else if(min_chi2<1000.0)
					chi[i-qini] = (intensity_exp[exp_index+(i-qini)]-diff_intensity[i])/pow(sigma[exp_index+(i-qini)]/sqrt(10.0),2);
			else 
					chi[i-qini] = (intensity_exp[exp_index+(i-qini)]-diff_intensity[i])/pow(sigma[exp_index+(i-qini)],2);
			
			chi[i-qini] *= chi2_scale;
		
		}
	}

	chi2=0.0;

	for(i=qini ; i<=100 ; i++){
		chi2 += pow((intensity_exp[exp_index+(i-qini)]-diff_intensity[i]),2)/pow(sigma[exp_index+(i-qini)],2);
	}
	chi2 *= chi2_scale;
	
	if(indicator == true){
		out.open("./result/chi2_rotation.dat",ios::app);
		static int count_step = 1;

		out<<count_step++<<"step : "<<chi2<<endl;
		out.close();
	}

	else{
		out.open("./result/chi2_RBMD.dat", ios::app);
		static int count_step2 = 1;
		out<<count_step2++<<"step : "<<chi2<<endl;
		out.close();
	}
}

void MakeFile(){
	ofstream out;
	out.open("./result/para_result.dati",ios::app);

	out<<"parameter result"<<endl;for(int i=0 ; i<num_of_parameter ; i++){
		out<<"para_"<<i<<" : "<<para_result[i]<<endl;	
	}
	out.close();
	MakePDB();
}

//additional functions for RBMD *********************************************************

void InitialBondCheck(){
	cout<<"InitialBondCheck()"<<endl;

	d_initial = 0.0;
	double temp;
	int i;

	for(i=0 ; i<num_of_rb/2-1 ; i++){
		temp = (com_x[i]-x[HEM_A_NC])*(com_x[i]-x[HEM_A_NC]) + (com_y[i]-y[HEM_A_NC])*(com_y[i]-y[HEM_A_NC]) + (com_z[i]-z[HEM_A_NC])*(com_z[i]-z[HEM_A_NC]);
		temp = sqrt(temp);
		d_initial += temp;
	}

	for(i=num_of_rb/2 ; i<num_of_rb-1 ; i++){
		temp = (com_x[i]-x[HEM_B_NC])*(com_x[i]-x[HEM_B_NC]) + (com_y[i]-y[HEM_B_NC])*(com_y[i]-y[HEM_B_NC]) + (com_z[i]-z[HEM_B_NC])*(com_z[i]-z[HEM_B_NC]);
		temp = sqrt(temp);
		d_initial += temp;
    }
}

void ReadFormfactor(){
	cout<<"ReadFormfactor()"<<endl;
	int i,j;
	ifstream in;
	in.open(formfactor);
	if(in==NULL){
		cerr<<"formfactor file is not opened."<<endl;
		exit(1);
	}

	for(i=0 ; i<num_of_atom_protein ; i++){
		for(j=0 ; j<=Nq-qini ; j++){
			in>>form[i][j];
		}
	}
	in.close();
}

void ReadMass(){
	cout<<"ReadMass()"<<endl;
	int i;
	ifstream in;
	in.open(reducedmass);
	if(in==NULL){
		cerr<<"mass file is not opened."<<endl;
		exit(1);
	}

	for(i=0 ; i<num_of_atom_protein ; i++){
		in>>mass[i];
	}
	in.close();

	for(i=0 ; i<num_of_rb ; i++)
		mass_rb[i] = 0.0;
		
	int temp_rb_index=0;

	for(i=0 ; i<num_of_atom_protein ; i++){
		if(i<=fin_rb[temp_rb_index]){
			mass_rb[temp_rb_index] += mass[i];
			if(i==fin_rb[temp_rb_index])
				temp_rb_index++;
		}
	}

//	for(i=0 ; i<num_of_rb ; i++)
//		cout<<mass_rb[i]<<endl;
}

void CalCom(){
	int i;
	
	for(i=0 ; i<num_of_rb ; i++)
		com_x[i] = com_y[i] = com_z[i] = 0.0;

	int temp_rb_index = 0;

	for(i=0 ; i<num_of_atom_protein ; i++){
		if(i<=fin_rb[temp_rb_index]){
			com_x[temp_rb_index] += (x[i]*mass[i]);
			com_y[temp_rb_index] += (y[i]*mass[i]);
			com_z[temp_rb_index] += (z[i]*mass[i]);
		}
		if(i==fin_rb[temp_rb_index])
			temp_rb_index++;
	}

	for(i=0 ; i<num_of_rb ; i++){
		com_x[i] /= mass_rb[i];
		com_y[i] /= mass_rb[i];
		com_z[i] /= mass_rb[i];
	}
} // calculate COM of each RB

void CalRelativePosition(){
	int temp_rb_index = 0;
	int i;

	for(i=0 ; i<num_of_atom_protein ; i++){
		if(i<=fin_rb[temp_rb_index]){
			x_rp[i] = x[i] - com_x[temp_rb_index];
			y_rp[i] = y[i] - com_y[temp_rb_index];
			z_rp[i] = z[i] - com_z[temp_rb_index];
		}
		if(i==fin_rb[temp_rb_index])
			temp_rb_index++;
	}
}

void CalForce(){
	cout<<"CalForce()"<<endl;
	int i,j,k;

	double temp_r, temp_r2;
	double temp_x, temp_y, temp_z, temp_x2, temp_y2, temp_z2;
	double temp_force;
	double temp_fx, temp_fy, temp_fz;
	double q;
	int kkk=0;
	int times = 0;
	for(i=0 ; i<num_of_atom_protein ; i++)
		f_lj_x[i] = f_lj_y[i] = f_lj_z[i] = f_chi_x[i] = f_chi_y[i] = f_chi_z[i] = f_sym_x[i] = f_sym_y[i] = f_sym_z[i] = 0.0;
	
	bool in_loop = true; 
	int ix;

	for(k=0 ; k<num_of_rb-1 ; k++){
		for(i=ini_rb[k] ; i<=fin_rb[k] ; i++){
			for(j=ini_rb[k+1] ; j<=fin_rb[num_of_rb-1] ; j++){
				for(ix = 0 ; ix<num_of_loop ; ix++){
					in_loop = true;
					if((i>=ini_loop[ix])&&(i<=fin_loop[ix]))
						in_loop = false;
					if((j>=ini_loop[ix])&&(j<=fin_loop[ix]))
						in_loop = false;
				}
						
				if(in_loop){		
					temp_x = x[i] - x[j];
					temp_y = y[i] - y[j];
					temp_z = z[i] - z[j];
					temp_r = sqrt(pow(temp_x, 2.0)+pow(temp_y, 2.0)+pow(temp_z, 2.0));


					if(temp_r > 0){			
						for(int x=0 ; x<=100-qini ; x++){
							q = (qini+x)*0.01;
							temp_force = p2*2*chi[x]*form[i][x]*form[j][x]*(cos(q*temp_r)-((sin(q*temp_r))/(q*temp_r)));
							temp_fx = temp_force*temp_x/temp_r/temp_r;
							temp_fy = temp_force*temp_y/temp_r/temp_r;
							temp_fz = temp_force*temp_z/temp_r/temp_r;
							f_chi_x[i] += temp_fx;
							f_chi_y[i] += temp_fy;
							f_chi_z[i] += temp_fz;
							f_chi_x[j] -= temp_fx;
							f_chi_y[j] -= temp_fy;
							f_chi_z[j] -= temp_fz;
						} // calculate chi force
					}
				}
			}
		}
	} 

	for(k=0 ; k<17 ; k++){
		for(i=ini_rb[k] ; i<=fin_rb[k] ; i++){
			for(j=ini_rb[k+1] ; j<=fin_rb[num_of_rb-1] ; j++){
				temp_x = x[i] - x[j];
				temp_y = y[i] - y[j];
				temp_z = z[i] - z[j];
				temp_r = sqrt(pow(temp_x, 2.0)+pow(temp_y, 2.0)+pow(temp_z, 2.0));

				if((a[i][0]== 'C') && (a[j][0]== 'N'))
					sig = sigNC;
				else if((a[j][0]== 'C') && (a[i][0]== 'N')) 
					sig = sigNC;
				else if(!strcmp(a[i], "FE") || !strcmp(a[j], "FE"))
					sig = sigFE;
				else
					sig = sigA;

				if((temp_r > 0)&&(temp_r<sig)){
					
						temp_force = p1*((48*(pow(sig,12.0)/pow(temp_r,13.0))) - (24*(pow(sig, 6.0)/pow(temp_r, 7.0))));
						temp_fx = temp_force*temp_x/temp_r;
						temp_fy = temp_force*temp_y/temp_r;
						temp_fz = temp_force*temp_z/temp_r;
						f_lj_x[i] += temp_fx;
						f_lj_y[i] += temp_fy;
						f_lj_z[i] += temp_fz;
						f_lj_x[j] -= temp_fx;
						f_lj_y[j] -= temp_fy;
						f_lj_z[j] -= temp_fz;
					//	UE_lj += p1*4*((pow(sig, 12.0)/pow(temp_r, 12.0))-(pow(sig, 6.0)/pow(temp_r, 6.0)));
				}				
			}
		}
	} 

	for(k=0 ; k<num_of_rb/2-1 ; k++){
        for(i=ini_rb[k] ; i<=fin_rb[k] ; i++){
            for(j=ini_rb[k+1] ; j<=fin_rb[num_of_rb/2-1] ; j++){
                    temp_x = x[i] - x[j];
                    temp_y = y[i] - y[j];
                    temp_z = z[i] - z[j];
                    temp_r = sqrt(pow(temp_x, 2.0)+pow(temp_y, 2.0)+pow(temp_z, 2.0));

                    temp_x2 = x[i+subunit_atom] - x[j+subunit_atom];
                    temp_y2 = y[i+subunit_atom] - y[j+subunit_atom];
                    temp_z2 = z[i+subunit_atom] - z[j+subunit_atom];
                    temp_r2 = sqrt(pow(temp_x2, 2.0)+pow(temp_y2, 2.0)+pow(temp_z2, 2.0));

                    temp_force = p3*(temp_r2-temp_r);

                    f_sym_x[i] += temp_force*temp_x/temp_r;
                    f_sym_y[i] += temp_force*temp_y/temp_r;
                    f_sym_z[i] += temp_force*temp_z/temp_r;
                    f_sym_x[j] -= temp_force*temp_x/temp_r;
                    f_sym_y[j] -= temp_force*temp_y/temp_r;
                    f_sym_z[j] -= temp_force*temp_z/temp_r;

                    f_sym_x[i+subunit_atom] -= temp_force*temp_x2/temp_r2;
                    f_sym_y[i+subunit_atom] -= temp_force*temp_y2/temp_r2;
                    f_sym_z[i+subunit_atom] -= temp_force*temp_z2/temp_r2;
                    f_sym_x[j+subunit_atom] += temp_force*temp_x2/temp_r2;
                    f_sym_y[j+subunit_atom] += temp_force*temp_y2/temp_r2;
                    f_sym_z[j+subunit_atom] += temp_force*temp_z2/temp_r2;
            }
        }
     }
	
	for(i=0 ; i<num_of_atom_protein ; i++){
		f_tot_x[i] = f_tot_y[i] = f_tot_z[i] = 0.0;
		f_tot_x[i] = f_lj_x[i] + f_chi_x[i] + f_sym_x[i];
		f_tot_y[i] = f_lj_y[i] + f_chi_y[i] + f_sym_y[i];
		f_tot_z[i] = f_lj_z[i] + f_chi_z[i] + f_sym_z[i];
	}

	double temp_chi_tot_x, temp_chi_tot_y, temp_chi_tot_z;

	for(i=0 ; i<num_of_rb ; i++){
		f_rb_tot_x[i] = f_rb_tot_y[i] = f_rb_tot_z[i] = 0.0;
//		temp_chi_tot_x = temp_chi_tot_y = temp_chi_tot_z = 0.0;
	
		for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
			f_rb_tot_x[i] += f_tot_x[j];
			f_rb_tot_y[i] += f_tot_y[j];
			f_rb_tot_z[i] += f_tot_z[j];
//			temp_chi_tot_x += f_chi_x[j];
//			temp_chi_tot_y += f_chi_y[j];
//			temp_chi_tot_z += f_chi_z[j];
		}
//		cout<<"chi_tot_x : "<<temp_chi_tot_x<<"   chi_tot_y : "<<temp_chi_tot_y<<"   chi_tot_z : "<<temp_chi_tot_z<<endl;
	}
}

void Translation(){
	cout<<"Translation()"<<endl;
	static int it_v=0;
	int i;

	for(i=0 ; i<num_of_rb ; i++){
		v_x[i] += f_rb_tot_x[i]*dt/mass_rb[i];
		v_y[i] += f_rb_tot_y[i]*dt/mass_rb[i];
		v_z[i] += f_rb_tot_z[i]*dt/mass_rb[i];
		v_x[i] = v_x[i]*sqrt(temperature/ini_temperature);
		v_y[i] = v_y[i]*sqrt(temperature/ini_temperature);
		v_z[i] = v_z[i]*sqrt(temperature/ini_temperature);
	}



}

void Rotation_Mb(){
	cout<<"Rotation()"<<endl;
	int i,j,k;
	static int it_w=0;
	double det; // determinent for moment of inertia 
	double dwx, dwy, dwz; // angular velocity change
	double w; // temporary memory to store length of angular velocity of each RB
	double rduw; // temporary memory to store r_p dot unitvector of w
	double auwcr; // temporary memory to store length of (unit vector of w cross r_p)
	double rp2; // temporary memory to store square of r_p
	
	for(i=0 ; i<num_of_rb ; i++){
		torque_x[i] = torque_y[i] = torque_z[i] = 0.0;
		for(j=0 ; j<3 ; j++){
			for(k=0 ; k<3 ; k++){
				moment_of_inertia[i][j][k] = 0.0;
			}
		}
	}

	for(i=0 ; i<num_of_rb ; i++){
		for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
			torque_x[i] += (((y_rp[j]*f_tot_z[j]) - (z_rp[j]*f_tot_y[j]))*p);
			torque_y[i] += (((z_rp[j]*f_tot_x[j]) - (x_rp[j]*f_tot_z[j]))*p);
			torque_z[i] += (((x_rp[j]*f_tot_y[j]) - (y_rp[j]*f_tot_x[j]))*p);
			moment_of_inertia[i][0][0] += (mass[j] * (pow(y_rp[j], 2.0) + pow(z_rp[j], 2.0)));
			moment_of_inertia[i][0][1] -= (mass[j] * x_rp[j] * y_rp[j]);
			moment_of_inertia[i][0][2] -= (mass[j] * x_rp[j] * z_rp[j]);
			moment_of_inertia[i][1][0] -= (mass[j] * x_rp[j] * y_rp[j]);
			moment_of_inertia[i][1][1] += (mass[j] * (pow(x_rp[j], 2.0) + pow(z_rp[j], 2.0)));
			moment_of_inertia[i][1][2] -= (mass[j] * y_rp[j] * z_rp[j]);
			moment_of_inertia[i][2][0] -= (mass[j] * x_rp[j] * z_rp[j]);
			moment_of_inertia[i][2][1] -= (mass[j] * y_rp[j] * z_rp[j]);
			moment_of_inertia[i][2][2] += (mass[j] * (pow(x_rp[j], 2.0) + pow(y_rp[j], 2.0)));
		}
	}	// torque calculation for each RB and calculate moment of inertia for each RB


	for(i=0 ; i<num_of_rb ; i++){
		Iwx[i] = (moment_of_inertia[i][0][0]*w_x[i]) + (moment_of_inertia[i][0][1]*w_y[i]) + (moment_of_inertia[i][0][2]*w_z[i]);
		Iwy[i] = (moment_of_inertia[i][1][0]*w_x[i]) + (moment_of_inertia[i][1][1]*w_y[i]) + (moment_of_inertia[i][1][2]*w_z[i]);
		Iwz[i] = (moment_of_inertia[i][2][0]*w_x[i]) + (moment_of_inertia[i][2][1]*w_y[i]) + (moment_of_inertia[i][2][2]*w_z[i]);
		Ax[i] = torque_x[i] - ((w_y[i]*Iwz[i])-(w_z[i]*Iwy[i]));
		Ay[i] = torque_y[i] - ((w_z[i]*Iwx[i])-(w_x[i]*Iwz[i]));
		Az[i] = torque_z[i] - ((w_x[i]*Iwy[i])-(w_y[i]*Iwx[i]));
	

		det = moment_of_inertia[i][0][0]*moment_of_inertia[i][1][1]*moment_of_inertia[i][2][2] - moment_of_inertia[i][0][0]*moment_of_inertia[i][1][2]*moment_of_inertia[i][2][1] - moment_of_inertia[i][0][1]*moment_of_inertia[i][1][0]*moment_of_inertia[i][2][2] + moment_of_inertia[i][0][1]*moment_of_inertia[i][2][0]*moment_of_inertia[i][1][2] + moment_of_inertia[i][0][2]*moment_of_inertia[i][1][0]*moment_of_inertia[i][2][1] - moment_of_inertia[i][0][2]*moment_of_inertia[i][2][0]*moment_of_inertia[i][1][1];

		
		if(det!=0){
			moment_inverse[i][0][0] = (1/det)*((moment_of_inertia[i][1][1]*moment_of_inertia[i][2][2]) - (moment_of_inertia[i][1][2]*moment_of_inertia[i][2][1]));
			moment_inverse[i][0][1] = (1/det)*((moment_of_inertia[i][0][2]*moment_of_inertia[i][2][1]) - (moment_of_inertia[i][0][1]*moment_of_inertia[i][2][2]));
			moment_inverse[i][0][2] = (1/det)*((moment_of_inertia[i][0][1]*moment_of_inertia[i][1][2]) - (moment_of_inertia[i][0][2]*moment_of_inertia[i][1][1]));
			moment_inverse[i][1][0] = (1/det)*((moment_of_inertia[i][1][2]*moment_of_inertia[i][2][0]) - (moment_of_inertia[i][1][0]*moment_of_inertia[i][2][2]));
			moment_inverse[i][1][1] = (1/det)*((moment_of_inertia[i][0][0]*moment_of_inertia[i][2][2]) - (moment_of_inertia[i][0][2]*moment_of_inertia[i][2][0]));
			moment_inverse[i][1][2] = (1/det)*((moment_of_inertia[i][0][2]*moment_of_inertia[i][1][0]) - (moment_of_inertia[i][0][0]*moment_of_inertia[i][1][2]));
			moment_inverse[i][2][0] = (1/det)*((moment_of_inertia[i][1][0]*moment_of_inertia[i][2][1]) - (moment_of_inertia[i][1][1]*moment_of_inertia[i][2][0]));
			moment_inverse[i][2][1] = (1/det)*((moment_of_inertia[i][0][1]*moment_of_inertia[i][2][0]) - (moment_of_inertia[i][0][0]*moment_of_inertia[i][2][1]));
			moment_inverse[i][2][2] = (1/det)*((moment_of_inertia[i][0][0]*moment_of_inertia[i][1][1]) - (moment_of_inertia[i][0][1]*moment_of_inertia[i][1][0]));
		}
		
		dwx = (moment_inverse[i][0][0]*Ax[i] + moment_inverse[i][0][1]*Ay[i] + moment_inverse[i][0][2]*Az[i])*dt;
		dwy = (moment_inverse[i][1][0]*Ax[i] + moment_inverse[i][1][1]*Ay[i] + moment_inverse[i][1][2]*Az[i])*dt;
		dwz = (moment_inverse[i][2][0]*Ax[i] + moment_inverse[i][2][1]*Ay[i] + moment_inverse[i][2][2]*Az[i])*dt;
	
		w_x[i] += dwx;
		w_y[i] += dwy;
		w_z[i] += dwz;
		
		w_x[i] = w_x[i]*sqrt(temperature/ini_temperature);
		w_y[i] = w_y[i]*sqrt(temperature/ini_temperature);
		w_z[i] = w_z[i]*sqrt(temperature/ini_temperature); // temperature dependency of random velocity

	}

	for(i=0 ; i<num_of_rb ; i++){
		w = sqrt(pow(w_x[i], 2.0)+pow(w_y[i], 2.0)+pow(w_z[i], 2.0));

		if(w!=0){
			unit_w_x[i] = w_x[i]/w;
			unit_w_y[i] = w_y[i]/w;
			unit_w_z[i] = w_z[i]/w;
		} // unit vector of angular velocity 

		for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
			rduw = x_rp[j]*unit_w_x[i] + y_rp[j]*unit_w_y[i] + z_rp[j]*unit_w_z[i];
			uwcrx[j] = unit_w_y[i]*z_rp[j] - unit_w_z[i]*y_rp[j];
			uwcry[j] = unit_w_z[i]*x_rp[j] - unit_w_x[i]*z_rp[j];
			uwcrz[j] = unit_w_x[i]*y_rp[j] - unit_w_y[i]*x_rp[j];

			auwcr = sqrt(pow(uwcrx[j], 2.0)+pow(uwcry[j], 2.0)+pow(uwcrz[j], 2.0));
			rp2 =  pow(x_rp[j],2.0)+pow(y_rp[j],2.0)+pow(z_rp[j],2.0);

			if(auwcr!=0){
				x_rp[j] = rduw * unit_w_x[i] + cos(w*dt)*(x_rp[j] - (rduw * unit_w_x[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcrx[j];
				y_rp[j] = rduw * unit_w_y[i] + cos(w*dt)*(y_rp[j] - (rduw * unit_w_y[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcry[j];
				z_rp[j] = rduw * unit_w_z[i] + cos(w*dt)*(z_rp[j] - (rduw * unit_w_z[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcrz[j];
			}
		/*	if(i==0)
				cout<<"x_rp : "<<x_rp[j]<<endl;*/
		}
	}
}

void PositionUpdate(){
	cout<<"PositionUpdate()"<<endl;
        double w; // temporary memory to store length of angular velocity of each RB
        double rduw; // temporary memory to store r_p dot unitvector of w
        double auwcr; // temporary memory to store length of (unit vector of w cross r_p)
        double rp2; // temporary memory to store square of r_p

	int i,j;
	
	for(i=0 ; i<num_of_rb ; i++){
		com_x[i] += (v_x[i] * dt);
		com_y[i] += (v_y[i] * dt);
		com_z[i] += (v_z[i] * dt);
	}
	// calculate position of COM of each RB by translational velocity

	if(indicator_PositionUpdate=true){
		for(i=0 ; i<num_of_rb ; i++){
			if(w_x[i] == 0) goto here;

                w = sqrt(pow(w_x[i], 2.0)+pow(w_y[i], 2.0)+pow(w_z[i], 2.0));

                        unit_w_x[i] = w_x[i]/w;
                        unit_w_y[i] = w_y[i]/w;
                        unit_w_z[i] = w_z[i]/w;

                for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
                        rduw = x_rp[j]*unit_w_x[i] + y_rp[j]*unit_w_y[i] + z_rp[j]*unit_w_z[i];
                        uwcrx[j] = unit_w_y[i]*z_rp[j] - unit_w_z[i]*y_rp[j];
                        uwcry[j] = unit_w_z[i]*x_rp[j] - unit_w_x[i]*z_rp[j];
                        uwcrz[j] = unit_w_x[i]*y_rp[j] - unit_w_y[i]*x_rp[j];

                        auwcr = sqrt(pow(uwcrx[j], 2.0)+pow(uwcry[j], 2.0)+pow(uwcrz[j], 2.0));
                        rp2 =  pow(x_rp[j],2.0)+pow(y_rp[j],2.0)+pow(z_rp[j],2.0);

                        if(auwcr!=0){
                                x_rp[j] = rduw * unit_w_x[i] + cos(w*dt)*(x_rp[j] - (rduw * unit_w_x[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcrx[j];
                                y_rp[j] = rduw * unit_w_y[i] + cos(w*dt)*(y_rp[j] - (rduw * unit_w_y[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcry[j];
                                z_rp[j] = rduw * unit_w_z[i] + cos(w*dt)*(z_rp[j] - (rduw * unit_w_z[i])) + sin(w*dt)*(sqrt(rp2-pow(rduw,2))/auwcr)*uwcrz[j];
                        }
                /*      if(i==0)
 *                                      cout<<"x_rp : "<<x_rp[j]<<endl;*/
                }
			
        }	
	}


// remember relative position(before updated) : ( loop last 'C' - the next helix first atom ) ********************************
	double x_rel[num_of_loop], y_rel[num_of_loop], z_rel[num_of_loop];
	int index;
	
	for(i=0 ; i<num_of_loop/2 ; i++){
		index = 0;

		while(true){
			if(!strcmp(a[fin_loop[i]-index], "C")){
				x_rel[i] = x[fin_loop[i]-index];
				y_rel[i] = y[fin_loop[i]-index];
				z_rel[i] = z[fin_loop[i]-index];
				break;
			}
			index++;
		}

		x_rel[i] -= x[ini_rb[i+1]];
		y_rel[i] -= y[ini_rb[i+1]];
		z_rel[i] -= z[ini_rb[i+1]];
	}


	for(i=num_of_loop/2 ; i<num_of_loop ; i++){
		index = 0;

		while(true){
			if(!strcmp(a[fin_loop[i]-index], "C")){
				x_rel[i] = x[fin_loop[i]-index];
				y_rel[i] = y[fin_loop[i]-index];
				z_rel[i] = z[fin_loop[i]-index];
				break;
			}       
			index++;
		}

		x_rel[i] -= x[ini_rb[i+3]];
		y_rel[i] -= y[ini_rb[i+3]];
		z_rel[i] -= z[ini_rb[i+3]];
        } 

//position update (loop & helix) ***********************************************************
here:
	for(i=0 ; i<num_of_rb ; i++){
		for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
			x[j] = com_x[i] + x_rp[j];
			y[j] = com_y[i] + y_rp[j];
			z[j] = com_z[i] + z_rp[j];
		
		/*if(i==0)
			cout<<j<<" : updated position x : "<<x[j]<<endl;*/
		}
	}

// His101 - Fe rearrangement ***************************************************************
	if(w_x[0]!=0){
		double his_x, his_y, his_z;
		double len;
		double dot1, dot2;
				
		his_x = x[Fe_A]-x[HIS_101_NE2_A];
		his_y = y[Fe_A]-y[HIS_101_NE2_A];
		his_z = z[Fe_A]-z[HIS_101_NE2_A];

		len = sqrt(pow(his_x,2)+pow(his_y,2)+pow(his_z,2));
		
		if(len > Fe_HIS_101){
			dot1 = (f_rb_tot_x[HIS_101_rb]*his_x) + (f_rb_tot_y[HIS_101_rb]*his_y) + (f_rb_tot_z[HIS_101_rb]*his_z);
			dot2 = (f_rb_tot_x[Fe_rb]*his_x) + (f_rb_tot_y[Fe_rb]*his_y) + (f_rb_tot_z[Fe_rb]*his_z);
			
			double temp_len =  len - Fe_HIS_101;
			if((dot1!=0)&&(dot2!=0)){
			com_x[HIS_101_rb] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_x/len;
			com_y[HIS_101_rb] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_y/len;
			com_z[HIS_101_rb] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_z/len;

			com_x[Fe_rb] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_x/len;
			com_y[Fe_rb] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_y/len;
			com_z[Fe_rb] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_z/len;
			}
		}

		his_x = x[Fe_B]-x[HIS_101_NE2_B];
		his_y = y[Fe_B]-y[HIS_101_NE2_B];
		his_z = z[Fe_B]-z[HIS_101_NE2_B];

		len = sqrt(pow(his_x,2)+pow(his_y,2)+pow(his_z,2));
		
		
		
		if(len>Fe_HIS_101){
			dot1 = (f_rb_tot_x[HIS_101_rb_B]*his_x) + (f_rb_tot_y[HIS_101_rb_B]*his_y) + (f_rb_tot_z[HIS_101_rb_B]*his_z);
			dot2 = (f_rb_tot_x[Fe_rb_B]*his_x) + (f_rb_tot_y[Fe_rb_B]*his_y) + (f_rb_tot_z[Fe_rb_B]*his_z);
			
			double temp_len =  len - Fe_HIS_101;
			
			if((dot1!=0)&&(dot2!=0)){
			com_x[HIS_101_rb_B] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_x/len;
			com_y[HIS_101_rb_B] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_y/len;
			com_z[HIS_101_rb_B] += (fabs(dot1)/(fabs(dot1)+fabs(dot2)))*temp_len*his_z/len;

			com_x[Fe_rb_B] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_x/len;
			com_y[Fe_rb_B] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_y/len;
			com_z[Fe_rb_B] -= (fabs(dot2)/(fabs(dot1)+fabs(dot2)))*temp_len*his_z/len;
			}
		}

		for(i=0 ; i<num_of_rb ; i++){
			if((i==HIS_101_rb)||(i==Fe_rb)||(i==HIS_101_rb_B)||(i==Fe_rb_B)){
				for(j=ini_rb[i] ; j<=fin_rb[i] ; j++){
					x[j] = com_x[i] + x_rp[j];
					y[j] = com_y[i] + y_rp[j];
					z[j] = com_z[i] + z_rp[j];
				}
			}
		}
	}
// His101 - Fe rearrangement is finished. ************************************************** 


// Translation of amino acid ***************************************************************
	int residue;
	double x_loop, y_loop, z_loop; // the extent of loop movement from before updated to updated
	
	if(w_x[0]!=0){
		for(i=0 ; i<num_of_loop/2 ; i++){
			residue = 0;
			index = 0;
						
			x_loop = x[ini_rb[i+1]];
			y_loop = y[ini_rb[i+1]];
			z_loop = z[ini_rb[i+1]];
			while(true){
				if(!strcmp(a[fin_loop[i]-index], "C")){
					x_loop -= x[fin_loop[i]-index];
					y_loop -= y[fin_loop[i]-index];
					z_loop -= z[fin_loop[i]-index];
					break;
				}
				index++;
			}

			x_loop += x_rel[i]; 
			y_loop += y_rel[i]; 
			z_loop += z_rel[i];
						
			x_loop /= (num_in_loop[i]+1);
			y_loop /= (num_in_loop[i]+1);
			z_loop /= (num_in_loop[i]+1);		


			for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){

				if(!strcmp(a[j],"N"))
					residue++;

				x[j] += x_loop*residue;
				y[j] += y_loop*residue;
				z[j] += z_loop*residue;
			}
		}

		for(i=num_of_loop/2 ; i<num_of_loop ; i++){
			residue = 0;
			index = 0;
			x_loop = x[ini_rb[i+3]];
			y_loop = y[ini_rb[i+3]];
			z_loop = z[ini_rb[i+3]];
			while(true){
				if(!strcmp(a[fin_loop[i]-index], "C")){
					x_loop -= x[fin_loop[i]-index];
					y_loop -= y[fin_loop[i]-index];
					z_loop -= z[fin_loop[i]-index];
					break;
				}
				index++;
			} 
			x_loop += x_rel[i];
			y_loop += y_rel[i];
			z_loop += z_rel[i];

			x_loop /= (num_in_loop[i]+1);
			y_loop /= (num_in_loop[i]+1);
			z_loop /= (num_in_loop[i]+1);

			for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){

				if(!strcmp(a[j],"N"))
					residue++;

				x[j] += x_loop*residue;
				y[j] += y_loop*residue;
				z[j] += z_loop*residue;
			}
		}
	} 

// Translation of amino acid to maintain distance between adjacent two amino acids
/*
	double *x_temp, *y_temp, *z_temp, *r_temp;
	double ave_r;
	double max, max_index;

	for(i=0 ; i<14 ; i++){
//		if((i==0)||(i==6)||(i==7)||(i==13)){
		index = 0;
		residue = 0;
		ave_r = 0;
		max_index = 0;

		x_temp = new double[num_in_loop[i]+2];
		y_temp = new double[num_in_loop[i]+2];
		z_temp = new double[num_in_loop[i]+2];
		r_temp = new double[num_in_loop[i]+1];

		while(true){
			if(!strcmp(a[ini_loop[i]-index],"CA")){
				x_temp[0] = x[ini_loop[i]-index];
				y_temp[0] = y[ini_loop[i]-index];
				z_temp[0] = z[ini_loop[i]-index];
				break;
			}
			index++;
		}

		for(j=ini_loop[i] ; residue<num_in_loop[i]+1 ; j++){
			if(!strcmp(a[j],"CA")){
				residue++;
				x_temp[residue] = x[j];
				y_temp[residue] = y[j];
				z_temp[residue] = z[j];
			}
		}
		
		for(j=0 ; j<num_in_loop[i]+1 ; j++){
			r_temp[j] = sqrt(pow(x_temp[j+1]-x_temp[j],2)+pow(y_temp[j+1]-y_temp[j],2)+pow(z_temp[j+1]-z_temp[j],2));
			ave_r += r_temp[j];
		}
		ave_r /= (num_in_loop[i]+1);
		if(ave_r>3.85){
			cout<<"ave_r : "<<ave_r<<endl;

		for(j=0 ; j<num_in_loop[i]+1 ; j++){
			if(j==0){
				max = r_temp[j];
				max_index = 0;
			}
			else{
				if(max<r_temp[j]){
					max = r_temp[j];
					max_index = j;
				}
			}
		}
		double temp, ratio, dist;

		if(max_index == num_in_loop[i]){
			residue = 0;
			for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
				if(!strcmp(a[j],"N")){
					residue++;
					temp = max - ave_r;
					temp /= (num_in_loop[i]+1);
					temp *= residue;
					dist = sqrt(pow(x_temp[residue+1]-x_temp[residue],2) + pow(y_temp[residue+1]-y_temp[residue],2) +  pow(z_temp[residue+1]-z_temp[residue],2));
					ratio = temp/dist;
				}
				x[j] += ((x_temp[residue+1]-x_temp[residue])*ratio);
				y[j] += ((y_temp[residue+1]-y_temp[residue])*ratio);
				z[j] += ((z_temp[residue+1]-z_temp[residue])*ratio);
			}
		} // if the last bond is elongated, rearrange positions from the last.
		
		else{
			residue = 0;
			for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
				if(!strcmp(a[j], "N")){
					residue++;
					if(residue>max_index){
						temp = max - ave_r;
						temp /= (num_in_loop[i]-max_index+1);
						temp *= (num_in_loop[i]-residue+1);
						dist = sqrt(pow(x_temp[residue-1]-x_temp[residue],2) + pow(y_temp[residue-1]-y_temp[residue],2) +  pow(z_temp[residue-1]-z_temp[residue],2));
						ratio = temp/dist;
					}
				}
				if(residue>max_index){
					x[j] += ((x_temp[residue-1]-x_temp[residue])*ratio);
					y[j] += ((y_temp[residue-1]-y_temp[residue])*ratio);
					z[j] += ((z_temp[residue-1]-z_temp[residue])*ratio);
				}
			}
		}
		}
		delete x_temp;
		delete y_temp;
		delete z_temp;
		delete r_temp;
//	}
	}
*/
// If the bond is broken due to large loop change, then adjust orientation of amino acids.

	bool bond_check, bond_check_peptide;
	bond_check = BondBreak();
	bond_check_peptide = BondBreakPeptide();
	
//	for(int k=0 ; k<num_of_break ; k++){
//		cout<<num_break_rb[k]<<endl;
//	}

	//if(w_x[0]!=0)
	{
			if(bond_check==true){
				if(termi_break==true){
					BringTerminus();
					bond_check_peptide = BondBreakPeptide();
					if(bond_check_peptide==true)
						AArearrange();
				}

				if(inter_break==true){
					Looprearrange();
					bond_check_peptide = BondBreakPeptide();
					if(bond_check_peptide==true)
						AArearrange();
					bond_check = BondBreak();
					//			if(bond_check==true)
					//				ConstrainMove();
				}
			}			

			else if((bond_check==false)&&(bond_check_peptide==true))
				AArearrange();	
	}

/*
	// protocol 2.
	
	         if(bond_check==true)
	                         ConstrainMove();
	                                 if((bond_check==false)&&(bond_check_peptide==true))
										AArearrange();	
*/
}

void ConstrainMove(){
	cout<<"ConstrainMove()"<<endl;
	system("cp ./prev.pdb new.pdb");
	ReadPDB();
	CalCom();
	CalRelativePosition();
	
	int temp;

	for(int i=0 ; i<num_of_break ; i++){
//		if((num_break_rb[i]!=0)&&(num_break_rb[i]!=6)&&(num_break_rb[i]!=9)&&(num_break_rb[i]!=15)){
			temp = num_break_rb[i];
			v_x[temp] = v_y[temp] = v_z[temp] = w_x[temp] = w_y[temp] = w_z[temp] = 0.0;
			temp += 1;
			v_x[temp] = v_y[temp] = v_z[temp] = w_x[temp] = w_y[temp] = w_z[temp] = 0.0;
//		}
	}
	
	if(indicator_PositionUpdate==true)
		PositionUpdate();

	else{
		indicator_PositionUpdate = !indicator_PositionUpdate;
		PositionUpdate();
		indicator_PositionUpdate = !indicator_PositionUpdate;
	}	
}

// Rearrangement of amino acids in loops
void Looprearrange(){
	double *x_temp, *y_temp, *z_temp, *r_temp;
	double ave_r;
	double max, max_index;
	int i,j,index,residue;

	for(i=0 ; i<num_of_loop ; i++){
//		if((i==0)||(i==6)||(i==7)||(i==13)){
		index = 0;
		residue = 0;
		ave_r = 0;
		max_index = 0;

		x_temp = new double[num_in_loop[i]+2];
		y_temp = new double[num_in_loop[i]+2];
		z_temp = new double[num_in_loop[i]+2];
		r_temp = new double[num_in_loop[i]+1];

		while(true){
			if(!strcmp(a[ini_loop[i]-index],"CA")){
				x_temp[0] = x[ini_loop[i]-index];
				y_temp[0] = y[ini_loop[i]-index];
				z_temp[0] = z[ini_loop[i]-index];
				break;
			}
			index++;
		}//'CA' of rb which is the closest to the front of loop

		for(j=ini_loop[i] ; residue<num_in_loop[i]+1 ; j++){
			if(!strcmp(a[j],"CA")){
				residue++;
				x_temp[residue] = x[j];
				y_temp[residue] = y[j];
				z_temp[residue] = z[j];
			}
		}//'CA's in loop and the closest to the end of loop
		
		for(j=0 ; j<num_in_loop[i]+1 ; j++){
			r_temp[j] = sqrt(pow(x_temp[j+1]-x_temp[j],2)+pow(y_temp[j+1]-y_temp[j],2)+pow(z_temp[j+1]-z_temp[j],2));
			ave_r += r_temp[j];
		}
		ave_r /= (num_in_loop[i]+1);
		
		if(ave_r>3.85){
//			cout<<"ave_r : "<<ave_r<<endl;

			for(j=0 ; j<num_in_loop[i]+1 ; j++){
				if(j==0){
					max = r_temp[j];
					max_index = 0;
				}
				else{
					if(max<r_temp[j]){
						max = r_temp[j];
						max_index = j;
					}
				}
			}
			double temp, ratio, dist;

			if(max_index == num_in_loop[i]){
				residue = 0;
				for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
					if(!strcmp(a[j],"N")){
						residue++;
						temp = max - ave_r;
						temp /= (num_in_loop[i]+1);
						temp *= residue;
						dist = sqrt(pow(x_temp[residue+1]-x_temp[residue],2) + pow(y_temp[residue+1]-y_temp[residue],2) +  pow(z_temp[residue+1]-z_temp[residue],2));
						ratio = temp/dist;
					}
					x[j] += ((x_temp[residue+1]-x_temp[residue])*ratio);
					y[j] += ((y_temp[residue+1]-y_temp[residue])*ratio);
					z[j] += ((z_temp[residue+1]-z_temp[residue])*ratio);
				}
			} // if the bond between the last 2 'CA's of loop is elongated, AAs rearrange positions from the last.
			
			else{
				residue = 0;
				for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
					if(!strcmp(a[j], "N")){
						residue++;
						if(residue>max_index){
							temp = max - ave_r;
							temp /= (num_in_loop[i]-max_index+1);
							temp *= (num_in_loop[i]-residue+1);
							dist = sqrt(pow(x_temp[residue-1]-x_temp[residue],2) + pow(y_temp[residue-1]-y_temp[residue],2) +  pow(z_temp[residue-1]-z_temp[residue],2));
							ratio = temp/dist;
						}
					}
					if(residue>max_index){
						x[j] += ((x_temp[residue-1]-x_temp[residue])*ratio);
						y[j] += ((y_temp[residue-1]-y_temp[residue])*ratio);
						z[j] += ((z_temp[residue-1]-z_temp[residue])*ratio);
					}
				}
			}
		}
		delete x_temp;
		delete y_temp;
		delete z_temp;
		delete r_temp;
//	}
	}
}

// Rotation of amino acid in loops
void AArearrange(){
	int i,j,k;
	double vec_a[3], vec_b[3], vec_nor[3];
	int index;
	int residue;
	double *x_temp, *y_temp, *z_temp;
	double **b_temp;
	double cos_theta, cos_phi, sin_theta, sin_phi, len_nor, delta, cos_delta, sin_delta, len_nor_pro, delta_a, delta_b;
	double temp_x, temp_y, temp_z;
	double trans_x, trans_y, trans_z;

	for(i=0 ; i<num_of_loop ; i++){
		index = 0;
		residue = 0;

		x_temp = new double[num_in_loop[i]+2];
		y_temp = new double[num_in_loop[i]+2];
		z_temp = new double[num_in_loop[i]+2];
		b_temp = new double*[num_in_loop[i]+1];

		for(j=0 ; j<num_in_loop[i]+1 ; j++)
			b_temp[j] = new double[3];

		while(true){
			if(!strcmp(a[ini_loop[i]-index],"CA")){
				x_temp[0] = x[ini_loop[i]-index];
				y_temp[0] = y[ini_loop[i]-index];
				z_temp[0] = z[ini_loop[i]-index];
				break;
			}
			index++;
		}

		for(j=ini_loop[i] ; residue<num_in_loop[i]+1 ; j++){
			if(!strcmp(a[j],"CA")){
				residue++;
				x_temp[residue] = x[j];
				y_temp[residue] = y[j];
				z_temp[residue] = z[j];
			}
		}

		for(j=0 ; j<num_in_loop[i]+1 ; j++){
			b_temp[j][0] = x_temp[j+1]-x_temp[j];
			b_temp[j][1] = y_temp[j+1]-y_temp[j];
			b_temp[j][2] = z_temp[j+1]-z_temp[j];
		}

		residue = -1;
		int ini_num=0, fin_num=0;
		bool chk = true;

		double *x_restore, *y_restore, *z_restore;

		for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
			if(j==fin_loop[i]){
				chk = true;
				fin_num = j;
				goto aa_rot;
			}

			if(!strcmp(a[j],"N")){
				if(chk){
					chk = false;
					ini_num = j;
				}

				else{
					chk = true;
					fin_num = j-1;
					j--;
					goto aa_rot;
				}

				residue++;

				vec_a[0] = x[j];
				vec_a[1] = y[j];
				vec_a[2] = z[j];
				
				for(k=0 ; k<3 ; k++){
					vec_b[k] = b_temp[residue][k] + b_temp[residue+1][k];
					vec_b[k] /= 2;
				}
			}

			else if(!strcmp(a[j],"C")){
				vec_a[0] = x[j]-vec_a[0];
				vec_a[1] = y[j]-vec_a[1];
				vec_a[2] = z[j]-vec_a[2];
			}

			else if(!strcmp(a[j],"CA")){
				trans_x = x[j];
				trans_y = y[j];
				trans_z = z[j];
			}

			if(false){
aa_rot:				x_restore = new double[fin_num-ini_num+1];
					y_restore = new double[fin_num-ini_num+1];
					z_restore = new double[fin_num-ini_num+1];

				static int wi=1;
				wi++;
				double temp1, temp2;
				vec_nor[1] = 1.0;
				vec_nor[0] = (vec_a[1]*vec_b[2]-vec_b[1]*vec_a[2])/(vec_a[2]*vec_b[0]-vec_b[2]*vec_a[0]);
				vec_nor[2] = -(vec_a[0]*vec_nor[0]+vec_a[1])/vec_a[2];

				len_nor = sqrt(pow(vec_nor[0],2)+pow(vec_nor[1],2)+pow(vec_nor[2],2));
				len_nor_pro = sqrt(pow(vec_nor[0],2)+pow(vec_nor[2],2));

				cos_theta = vec_nor[1]/len_nor;
				cos_phi = vec_nor[2]/len_nor_pro;

				sin_theta = sqrt(1-pow(cos_theta,2));
				sin_phi = sqrt(1-pow(cos_phi,2));

// sin_phi sign check ****************************************************
				temp_x = vec_nor[0];
				temp_y = vec_nor[1];
				temp_z = vec_nor[2];
				vec_nor[0] = temp_x*cos_phi - temp_z*sin_phi;
				vec_nor[1] = temp_x*sin_theta*sin_phi + temp_y*cos_theta + temp_z*sin_theta*cos_phi;
				vec_nor[2] = temp_x*cos_theta*sin_phi - temp_y*sin_theta + temp_z*cos_theta*cos_phi;
				
/*				ofstream out;
				out.open("normal.dat",ios::app);
				out<<vec_nor[0]<<'\t'<<vec_nor[1]<<'\t'<<vec_nor[2]<<endl;
				out.close();
*/
				temp1 = vec_nor[0];

				sin_phi *= -1;

				vec_nor[0] = temp_x*cos_phi - temp_z*sin_phi;
				vec_nor[1] = temp_x*sin_theta*sin_phi + temp_y*cos_theta + temp_z*sin_theta*cos_phi;
				vec_nor[2] = temp_x*cos_theta*sin_phi - temp_y*sin_theta + temp_z*cos_theta*cos_phi;
				
				temp2 = vec_nor[0];
/*
				out.open("normal.dat",ios::app);
				out<<vec_nor[0]<<'\t'<<vec_nor[1]<<'\t'<<vec_nor[2]<<endl;
				out.close();
*/
				if(fabs(temp1)>fabs(temp2));
				else sin_phi*=-1;
				
// *************************************************************************

				int res_num = 0;
				for(k=ini_num ; k<=fin_num ; k++){
						x_restore[res_num] = x[k];
						y_restore[res_num] = y[k];
						z_restore[res_num] = z[k];
						res_num++;
				}


				for(k=ini_num ; k<=fin_num ; k++){
					x[k] -= trans_x;
					y[k] -= trans_y;
					z[k] -= trans_z;
				}
				
				for(k=ini_num ; k<=fin_num ; k++){
					temp_x = x[k];
					temp_y = y[k];
					temp_z = z[k];

					x[k] = temp_x*cos_phi - temp_z*sin_phi;
					y[k] = temp_x*sin_theta*sin_phi + temp_y*cos_theta + temp_z*sin_theta*cos_phi;
					z[k] = temp_x*cos_theta*sin_phi - temp_y*sin_theta + temp_z*cos_theta*cos_phi;
				}

				temp_x = vec_a[0];
				temp_y = vec_a[1];
				temp_z = vec_a[2];
				vec_a[0] = temp_x*cos_phi - temp_z*sin_phi;
				vec_a[1] = temp_x*sin_theta*sin_phi + temp_y*cos_theta + temp_z*sin_theta*cos_phi;
				vec_a[2] = temp_x*cos_theta*sin_phi - temp_y*sin_theta + temp_z*cos_theta*cos_phi;

				temp_x = vec_b[0];
				temp_y = vec_b[1];
				temp_z = vec_b[2];
				vec_b[0] = temp_x*cos_phi - temp_z*sin_phi;
				vec_b[1] = temp_x*sin_theta*sin_phi + temp_y*cos_theta + temp_z*sin_theta*cos_phi;
				vec_b[2] = temp_x*cos_theta*sin_phi - temp_y*sin_theta + temp_z*cos_theta*cos_phi;
				
				delta_a = acos(vec_a[0]/sqrt(pow(vec_a[0],2)+pow(vec_a[1],2)+pow(vec_a[2],2)));
				delta_b = acos(vec_b[0]/sqrt(pow(vec_b[0],2)+pow(vec_b[1],2)+pow(vec_b[2],2)));

				if(vec_a[2]<0)delta_a*=-1;
				if(vec_b[2]<0)delta_b*=-1;

				delta = fabs(delta_a-delta_b);

				sin_delta = sin(delta);
				cos_delta = cos(delta);

// delta sign check ***********************************************
				temp_x = vec_a[0];
				temp_z = vec_a[2];
				vec_a[0] = temp_x*cos_delta + temp_z*sin_delta;
				vec_a[2] = temp_z*cos_delta - temp_x*sin_delta;

				temp1 = (vec_a[0]*vec_b[0]) + (vec_a[1]*vec_b[1]) + (vec_a[2]*vec_b[2]);

				sin_delta *= -1;
	
				vec_a[0] = temp_x*cos_delta + temp_z*sin_delta;
				vec_a[2] = temp_z*cos_delta - temp_x*sin_delta;
				
				temp2 = (vec_a[0]*vec_b[0]) + (vec_a[1]*vec_b[1]) + (vec_a[2]*vec_b[2]);

				if(temp1>temp2)sin_delta *= -1;
// ****************************************************************
				

				for(k=ini_num ; k<=fin_num ; k++){
					temp_x = x[k];
					temp_z = z[k];
					x[k] = temp_x*cos_delta + temp_z*sin_delta;
					z[k] = temp_z*cos_delta - temp_x*sin_delta;
				}

				for(k=ini_num ; k<=fin_num ; k++){
					temp_x = x[k];
					temp_y = y[k];
					temp_z = z[k];
					x[k] = temp_x*cos_phi + temp_y*sin_phi*sin_theta + temp_z*sin_phi*cos_theta;
					y[k] = temp_y*cos_theta - temp_z*sin_theta;
					z[k] = temp_y*cos_phi*sin_theta - temp_x*sin_phi + temp_z*cos_phi*cos_theta;
				}

				for(k=ini_num ; k<=fin_num ; k++){
					x[k] += trans_x;
					y[k] += trans_y;
					z[k] += trans_z;
				}

				bool check = false;
				ofstream out;
				out.open("./nan_check");

				for(k=ini_num ; k<=fin_num ; k++){
						out<<x[k]<<endl<<y[k]<<endl<<z[k]<<endl;
				}

				out.close();
				ifstream in;
				in.open("./nan_check");
				char nan[255];

				for(k=0 ; k<(3*(fin_num-ini_num+1)) ; k++){
						in.getline(nan,255,'\n');
						if(!strcmp(nan,"nan") || !strcmp(nan,"-nan")){
								check = true;
								break;
						}
				}
				in.close();


				res_num = 0;

				if(check){
						for(k=ini_num ; k<=fin_num ; k++){
										x[k] = x_restore[res_num];
										y[k] = y_restore[res_num];
										z[k] = z_restore[res_num];
										res_num++;
						}
				}

				delete x_restore;
				delete y_restore;
				delete z_restore;



			}
			
		}
		cout<<residue<<endl;

		delete x_temp; 
		delete y_temp;
		delete z_temp;
		delete []b_temp;
	}
}

double BondCheck(){
	cout<<"BondCheck()"<<endl;
	double temp_d = 0.0;
	int i;
	double temp;

	for(i=0 ; i<num_of_rb/2-1 ; i++){
                temp = (com_x[i]-x[HEM_A_NC])*(com_x[i]-x[HEM_A_NC]) + (com_y[i]-y[HEM_A_NC])*(com_y[i]-y[HEM_A_NC]) + (com_z[i]-z[HEM_A_NC])*(com_z[i]-z[HEM_A_NC]);
                temp = sqrt(temp);
                temp_d += temp;
    }
        
	for(i=num_of_rb/2 ; i<num_of_rb-1 ; i++){
			temp = (com_x[i]-x[HEM_B_NC])*(com_x[i]-x[HEM_B_NC]) + (com_y[i]-y[HEM_B_NC])*(com_y[i]-y[HEM_B_NC]) + (com_z[i]-z[HEM_B_NC])*(com_z[i]-z[HEM_B_NC]);
			temp = sqrt(temp);
			temp_d += temp;
	}

	return temp_d;
}

double TargetFunction(){
	cout<<"TargetFunction()"<<endl;

// chi2 calculation ***************************************************
	cout<<"Chi2 calculation"<<endl;

	int i, j, k;

	ofstream out;
	out.open("diff.int");
	for(i=0 ; i<=Nq ; i++){
		diff_intensity[i] = intensity_final[i] - intensity_initial[i];
		diff_intensity[i] *= sf;
		//intensity_theory[i] = diff_intensity[i];
		out<<q[i]<<'\t'<<diff_intensity[i]<<endl;
	}
	out.close(); // calculate scattering intensity difference between two model 
		
	double temp;

	ifstream in;
	in.open(exp_curve);
	if(in==NULL)
	{
		cerr<<"exp_curve file is not open"<<endl;
		exit(1);
	}
	for(i=0 ; i<Nq-1 ; i++){
		in>>temp>>intensity_exp[i]>>sigma[i];
		if(temp == q[qini])
			exp_index = i;			
	}
	in.close();
	
	chi2=0.0;

	for(i=qini ; i<=100 ; i++){
		chi2 += pow((intensity_exp[exp_index+(i-qini)]-diff_intensity[i]),2)/pow(sigma[exp_index+(i-qini)],2);
	}
	chi2 *= chi2_scale;
	
	static int count_step2 = 0;

	if(indicator == true){
		out.open("./result/chi2_rotation.dat",ios::app);
		static int count_step = 1;

		out<<count_step++<<"step : "<<chi2<<endl;
		out.close();
	}

	else{
		out.open("./result/chi2_new_MD.dat", ios::app);
		out<<count_step2++<<"step : "<<chi2<<endl;
		out.close();
	}

	out.open("./result/target_function_component.dat", ios::app);
	sprintf(chi2_corrected, "%06.4f", chi2);
	out<<count_step2<<"step : chi2 : "<<chi2_corrected<<'\t';
	out.close();

//*********************************************************************************

// Penalty 1. retain distance between two RBs*******************************************

	penalty_RB_dist = 0.0;
	int index = 0;
	double temp_r;

	for(i=1 ; i<num_of_rb ; i++){
		while(true){
			if(!strcmp(a[fin_rb[i-1]-index],"CA")){
				temp_r = pow((x[ini_rb[i]+1]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]+1]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]+1]-z[fin_rb[i-1]-index]),2);
				temp_r = sqrt(temp_r);
				break;
			}
			index++;
		}
		index = 0;
		if(temp_r > 3.8){
			penalty_RB_dist += log(temp_r/3.8);
		}
	}	
	penalty_RB_dist *= fp1;
	
	out.open("./result/target_function_component.dat", ios::app);
	out<<"penalty_RB_dist : "<<penalty_RB_dist<<'\t';
	out.close();

//************************************************************************************

// Penalty 2. avoid collision*********************************************************
	penalty_collision = 0.0;

	double temp_r2;
	double temp_x, temp_y, temp_z, temp_x2, temp_y2, temp_z2;
	double temp_force;
	double temp_fx, temp_fy, temp_fz;
	double q;

	
	for(k=0 ; k<num_of_rb-1 ; k++){
		for(i=ini_rb[k] ; i<=fin_rb[k] ; i++){
			for(j=ini_rb[k+1] ; j<=fin_rb[num_of_rb-1] ; j++){
				temp_x = x[i] - x[j];
				temp_y = y[i] - y[j];
				temp_z = z[i] - z[j];
				temp_r = sqrt(pow(temp_x, 2.0)+pow(temp_y, 2.0)+pow(temp_z, 2.0));

				if((a[i][0]== 'C') && (a[j][0]== 'N'))
					sig = sigNC;
				else if((a[j][0]== 'C') && (a[i][0]== 'N')) 
					sig = sigNC;
				else if(!strcmp(a[i], "FE") || !strcmp(a[j], "FE"))
					sig = sigFE;
				else
					sig = sigA;

				if(temp_r > 0){
					if(temp_r < sig){
						penalty_collision += ((pow(sig,12.0)/pow(temp_r,12.0)) - (pow(sig, 6.0)/pow(temp_r, 6.0)));						
					} 
				}				
			}
		}
	} 

	penalty_collision *= fp2;

	out.open("./result/target_function_component.dat", ios::app);
	out<<"penalty_collision : "<<penalty_collision<<'\t';
	out.close();
//*********************************************************************************************************************

// Penalty 3. Symmetry of two subunit**********************************************************************************
	penalty_symmetry = 0.0;

	for(k=0 ; k<num_of_rb/2-1 ; k++){
		for(i=ini_rb[k] ; i<=fin_rb[k] ; i++){
			for(j=ini_rb[k+1] ; j<=fin_rb[num_of_rb/2-1] ; j++){
				temp_x = x[i] - x[j];
				temp_y = y[i] - y[j];
				temp_z = z[i] - z[j];
				temp_r = sqrt(pow(temp_x, 2.0)+pow(temp_y, 2.0)+pow(temp_z, 2.0));

				temp_x2 = x[i+subunit_atom] - x[j+subunit_atom];
				temp_y2 = y[i+subunit_atom] - y[j+subunit_atom];
				temp_z2 = z[i+subunit_atom] - z[j+subunit_atom];
				temp_r2 = sqrt(pow(temp_x2, 2.0)+pow(temp_y2, 2.0)+pow(temp_z2, 2.0));	
				
				penalty_symmetry += ((temp_r-temp_r2)*(temp_r-temp_r2));
			}
		}
	}

	penalty_symmetry *= fp3;

	out.open("./result/target_function_component.dat", ios::app);
	out<<"penalty_symmetry : "<<penalty_symmetry<<'\t';
	out.close();
//*********************************************************************************************************************

// Penalty 4. avoid unfolding******************************************************************************************
	d_final = BondCheck();

	penalty_unfolding = 0.0;

	if(d_final>d_initial)
		penalty_unfolding = ((d_final-d_initial)*(d_final-d_initial));

	penalty_unfolding *= fp4;

	out.open("./result/target_function_component.dat", ios::app);
	out<<"penalty_unfolding : "<<penalty_unfolding<<endl;
	out.close();
//*********************************************************************************************************************

	double temp_value = 0.0;

	temp_value = chi2 + penalty_factor*(penalty_RB_dist+penalty_collision+penalty_symmetry+penalty_unfolding);

	out.open("./result/target_value.dat", ios::app);
	out<<count_step2<<"step : target_value : "<<temp_value<<endl;
	out.close();

	return temp_value;
}

double CalRfactor(){
	int i;
	double R = 0.0;
	double R_up = 0.0;
	double R_down = 0.0;
		
	for(i=qini ; i<=100 ; i++){
		R_up += fabs(intensity_exp[exp_index+(i-qini)]-diff_intensity[i])/sigma[exp_index+(i-qini)];
	}

	for(i=qini ; i<=100 ; i++){
		R_down += fabs(intensity_exp[exp_index+(i-qini)])/sigma[exp_index+(i-qini)];
	}

	R = R_up/R_down;

	return R;
}

bool BondBreak(){
	bool bond_break = false;
	double temp_r;
	int index = 0;
	int i,j;

	for(i=0 ; i<num_of_rb ; i++)
		num_break_rb[i] = 0;

	num_of_break = 0;

	termi_break = false;
	inter_break = false;




	for(i=1 ; i<num_of_rb/2-1 ; i++){
		while(true){
			if(!strcmp(a[fin_rb[i-1]-index],"CA")){
				temp_r = pow((x[ini_rb[i]+1]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]+1]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]+1]-z[fin_rb[i-1]-index]),2);
				temp_r = sqrt(temp_r);
				break;
			}
			index++;
		}
		index = 0;
		if(temp_r > 4.6){
			num_break_rb[num_of_break++]=(i-1);
			if((i==1)||(i==num_of_rb/2-2))
				termi_break = true;
			else 
				inter_break = true;
			bond_break = true;
		}
	}

	index = 0;

	for(i=num_of_rb/2+1 ; i<num_of_rb-1 ; i++){
        while(true){
                if(!strcmp(a[fin_rb[i-1]-index],"CA")){
                        temp_r = pow((x[ini_rb[i]+1]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]+1]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]+1]-z[fin_rb[i-1]-index]),2);
                        temp_r = sqrt(temp_r);
                        break;
                }
                index++;
        }
        index = 0;
      		  if(temp_r > 4.6){
			num_break_rb[num_of_break++]=(i-1);
			if((i==num_of_rb/2+1)||(i==num_of_rb-2))
				termi_break = true;
			else
				inter_break = true;
        	        bond_break = true;
		}
        }

	double x1,y1,z1,x2,y2,z2;
	int temp_resi_num;
	bool temp_bool ;

	for(i=0 ; i<num_of_loop ; i++){
		index = 1;
		temp_bool = true;

		int pre_index=1;
		
		while(true){
			if(!strcmp(a[ini_loop[i]-pre_index],"CA"))
				break;		
			pre_index++;
		}

		for(j=ini_loop[i]-pre_index ; j<=fin_loop[i] ; j++){
			if(!strcmp(a[j],"CA")){
				if(index==1){
					x1 = x[j];
					y1 = y[j];
					z1 = z[j];
					index = 2;
				}
				else if(index==2){
					x2 = x[j];
					y2 = y[j];
					z2 = z[j];
					index = 3;
				}
			}

			if(index ==3){
				temp_r = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
				temp_r = sqrt(temp_r);

				if(temp_r>4.6){
					if(i>6) temp_resi_num = (i+2);
					else temp_resi_num = i;					
					
					for(int k=0 ; k<num_of_break ; k++){
						if(num_break_rb[k]==temp_resi_num)
							temp_bool=false;
					}
					
					if(temp_bool) num_break_rb[num_of_break++] = temp_resi_num;

					if((i==0)||(i==6)||(i==7)||(i==13))
						termi_break = true;
					else
						inter_break = true;
					bond_break = true;
				}
				x1 = x2;
				y1 = y2;
				z1 = z2;
				index=2;
			}
		}
	}

	return bond_break;
}

void BringTerminus(){
	cout<<"BringTeminus()"<<endl;
	int index = 0;
	int i,j,k;
	double x1,y1,z1,x2,y2,z2,temp_r;
	double direction[3];
	int index_break;
	int helper;
	double temp, len_direc;

	for(i=0 ; i<num_of_loop ; i++){
		if((i==0)||(i==6)||(i==7)||(i==13)){
			while(true){
				if(!strcmp(a[ini_loop[i]-index],"CA")){
					temp_r = pow((x[ini_loop[i]+1]-x[ini_loop[i]-index]),2) + pow((y[ini_loop[i]+1]-y[ini_loop[i]-index]),2) + pow((z[ini_loop[i]+1]-z[ini_loop[i]-index]),2);
					temp_r = sqrt(temp_r);
					break;	
				}
				index++;
			}

			if(temp_r > 4.6){
                                direction[0] = x[ini_loop[i]+1]-x[ini_loop[i]-index];
                                direction[1] = y[ini_loop[i]+1]-y[ini_loop[i]-index];
                                direction[2] = z[ini_loop[i]+1]-z[ini_loop[i]-index];

                                if((i==1)||(i==7))
                                        index_break = fin_rb[i-1]-index;

                                else if((i==6)||(i==13))
                                        index_break = ini_rb[i]+1;

                                if((i==6)||(i==13)){
                                        for( k=0 ; k<3 ; k++)
                                                direction[k] *= -1;
                                }

                                len_direc = sqrt(pow(direction[0],2)+pow(direction[1],2)+pow(direction[2],2));

				if(i==1){
					temp = temp_r - target_dist;
			
					for(k=ini_rb[0] ; k<=ini_loop[i]-1 ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
					}
                }
				
				else if(i==7){
					temp = temp_r - target_dist;
				
					for(k=ini_rb[9] ; k<=ini_loop[i]-1 ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
					}
				}

				else if(i==6){
					temp = temp_r - target_dist;

					for(k=ini_loop[i] ; k<=fin_rb[7] ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
					}
				}
				
				else if(i==13){
					temp = temp_r = target_dist;

					for(k=ini_loop[i] ; k<=fin_rb[16] ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
                    }
				}
            }
                        index = 0;
		}
	}


	index = 0;
	for(i=1 ; i<17 ; i++){
		if((i==1)||(i==7)||(i==10)||(i==16)){
			while(true){
				if(!strcmp(a[fin_rb[i-1]-index],"CA")){
					temp_r = pow((x[ini_rb[i]+1]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]+1]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]+1]-z[fin_rb[i-1]-index]),2);
					temp_r = sqrt(temp_r);
					break;
				}
				index++;
			}

			if(temp_r > 4.6){
				direction[0] = x[ini_rb[i]+1]-x[fin_rb[i-1]-index];
				direction[1] = y[ini_rb[i]+1]-y[fin_rb[i-1]-index];
				direction[2] = z[ini_rb[i]+1]-z[fin_rb[i-1]-index];
				
				if((i==1)||(i==10))
					index_break = fin_rb[i-1]-index;

				else if((i==7)||(i==16))
					index_break = ini_rb[i]+1;
				
				if((i==7)||(i==16)){
					for( k=0 ; k<3 ; k++)
						direction[k] *= -1;
				}

				len_direc = sqrt(pow(direction[0],2)+pow(direction[1],2)+pow(direction[2],2));

				helper = 0;

				if((i==1)||(i==10)){
					while(true){
						if(!strcmp(a[index_break+helper],"N"))
							break;
						helper++;;
					}

					index_break = index_break + helper-1;
					temp = temp_r - target_dist;
					
					for(k=ini_rb[i-1] ; k<=index_break ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
					}
				}

				else if((i==7)||(i==16)){
					index_break -= 1;					
					temp = temp_r - target_dist;

					for(k=index_break ; k<=fin_rb[i] ; k++){
						x[k] += temp*direction[0]/len_direc;
						y[k] += temp*direction[1]/len_direc;
						z[k] += temp*direction[2]/len_direc;
					}
				}		
			}
			index = 0;
		}
	}

	index = 0;

	for(i=0 ; i<num_of_loop ; i++){
		if((i==0)||(i==6)||(i==7)||(i==13)){
			index = 1;
			for(j=ini_loop[i] ; j<=fin_loop[i] ; j++){
				if(!strcmp(a[j],"CA")){
					if(index==1){
						x1 = x[j];
						y1 = y[j];
						z1 = z[j];
						index = 2;
					}
					else if(index==2){
						x2 = x[j];
						y2 = y[j];
						z2 = z[j];
						index = 3;
					}
				}

				if(index ==3){
					temp_r = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
					temp_r = sqrt(temp_r);

					if(temp_r>4.6){
						direction[0] = x2-x1;
						direction[1] = y2-y1;
						direction[2] = z2-z1;

						if((i==0)||(i==7))
							index_break = j;

						else if((i==6)||(i==13))
							index_break = j;
						
						if((i==6)||(i==13)){
							for( k=0 ; k<3 ; k++)
								direction[k] *= -1;
						}

						len_direc = sqrt(pow(direction[0],2)+pow(direction[1],2)+pow(direction[2],2));

						helper = 0;

						if((i==0)||(i==7)){
							index_break -= 2;
							temp = temp_r - target_dist;
							
							if(i==0){
								for(k=ini_rb[0] ; k<=index_break ; k++){
									x[k] += temp*direction[0]/len_direc;
									y[k] += temp*direction[1]/len_direc;
									z[k] += temp*direction[2]/len_direc;
								}
							}

							else if(i==7){
								for(k=ini_rb[9] ; k<=index_break ; k++){
									x[k] += temp*direction[0]/len_direc;
									y[k] += temp*direction[1]/len_direc;
									z[k] += temp*direction[2]/len_direc;
								}
							}
						}

						else if((i==6)||(i==13)){
							index_break -= 1;					
							temp = temp_r - target_dist;

							if(i==6){
								for(k=index_break ; k<=fin_rb[7] ; k++){
									x[k] += temp*direction[0]/len_direc;
									y[k] += temp*direction[1]/len_direc;
									z[k] += temp*direction[2]/len_direc;
								}
							}

							else if(i==13){
								for(k=index_break ; k<=fin_rb[16] ; k++){
									x[k] += temp*direction[0]/len_direc;
									y[k] += temp*direction[1]/len_direc;
									z[k] += temp*direction[2]/len_direc;
								}
							}
						}	
					}
					x1 = x2;
					y1 = y2;
					z1 = z2;
					index=2;
				}
			}
		}
	}
	cout<<"BringTerminus() finished"<<endl;
}

bool BondBreakPeptide(){
	bool bond_break = false;
	double temp_r;
	int index = 0;
	int i,j;
	double limit = 1.95;

	for(i=1 ; i<8 ; i++){
		while(true){
			if(!strcmp(a[fin_rb[i-1]-index],"C")){
				temp_r = pow((x[ini_rb[i]]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]]-z[fin_rb[i-1]-index]),2);
				temp_r = sqrt(temp_r);
				break;
			}
			index++;
		}
		index = 0;
		if(temp_r > limit)
			bond_break = true;
	}

	index = 0;

	for(i=10 ; i<17 ; i++){
        while(true){
                if(!strcmp(a[fin_rb[i-1]-index],"C")){
                        temp_r = pow((x[ini_rb[i]]-x[fin_rb[i-1]-index]),2) + pow((y[ini_rb[i]]-y[fin_rb[i-1]-index]),2) + pow((z[ini_rb[i]]-z[fin_rb[i-1]-index]),2);
                        temp_r = sqrt(temp_r);
                        break;
                }
                index++;
        }
        index = 0;
        if(temp_r > limit)
                bond_break = true;
    }

	double x1,y1,z1,x2,y2,z2;
	index = 0;

	for(i=0 ; i<num_of_loop ; i++){
		int pre_index=1;

                while(true){
                        if(!strcmp(a[ini_loop[i]-pre_index],"C"))
                                break;
                        pre_index++;
                }

		for(j=ini_loop[i]-pre_index ; j<=fin_loop[i]+1 ; j++){
			if(!strcmp(a[j],"C")){					
				x1 = x[j];
				y1 = y[j];
				z1 = z[j];
			}

			else if(!strcmp(a[j],"N")){
				x2 = x[j];
				y2 = y[j];
				z2 = z[j];
				index = 1;
			}

			if(index ==1){
				temp_r = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
				temp_r = sqrt(temp_r);

				if(temp_r>limit)
					bond_break = true;

				index=0;
			}
		}
	}

	return bond_break;
}

void SelectReplica(){
	cout<<"SelectReplica()"<<endl;
	
	first_time = true;
	least_chi2_index = 0;
	min_chi2 = 1e10;
	
	for(int i=0 ; i<1000 ; i++){
		least_chi2_table[i] = 1e9;
	}

	minimum = 1e10;

        ofstream out;
        out.open("./result/target_value_real.dat",ios::app);
        out<<"Select Replica"<<endl;

        system("ls replica >list.txt");
        ifstream in;
        in.open("list.txt");
        int num_of_line = -1;

        while(!in.eof()){
                in.getline(replica_line,500,'\n');
                num_of_line++;
        }

		if(num_of_line == 0)
		{
			cout << "program ended" << endl;
			exit(1);
		}
	
        in.close();

        ifstream in2;
        in2.open("list.txt");

        int random_num;
        random_num = random(num_of_line);

        for(int i=0 ; i<random_num ; i++)
                in2.getline(replica_line,500,'\n');

        in2.getline(replica_line,500,'\n');
        out<<replica_line<<" is selected."<<endl;
		cout<<replica_line<<" is selected."<<endl;
		in2.close();
        char cmd[500] = "mv ./replica/";
        strcat(cmd,replica_line);
        strcat(cmd," ./prev.pdb");
        system(cmd);
		cout<<cmd;
        system("rm list.txt");
        out.close();
        system("cp prev.pdb ./new.pdb");
        ReadPDB();
        CalCom();
        CalRelativePosition();
}
