/* 
** routines.h
**
** Header file for routines.c
*/

void initialise_parameters(SI (*));
void initialise_black_hole(PARTICLE (*));
void check_main_parameters(SI (*));
void calculate_parameters(SI (*)); 
void initialise_gridr(PARTICLE (*), SI (*));
void initialise_griddf(SI (*));
void initialise_haloStructure(SI (*));
void initialise_starStructure(SI (*));
void set_haloPositions(SI (*)); 
void set_starPositions(SI (*));
void set_haloVelocities(SI (*));
void set_starVelocities(SI (*));
void double_haloParticles(SI (*));
void double_starParticles(SI (*));
void set_index(SI (*));
void set_attributes(SI (*));
void calculate_stuff(PARTICLE (*), SI (*));
void usage(void);
void transfer_particles(const PARTICLE (*),const SI (*), TIPSY_STRUCTURE (*));
void select_stars(SI (*));
void set_lighttomass(SI (*));
void displace(SI (*));
void rotate(void);
