/* 
** io.h
**
** Header file for io.h
*/

void write_griddf( FILE (*),const SI (*));
void write_gridr(FILE (*), const SI (*));
void write_ifrit(FILE (*),const SI (*));
void write_profiles(FILE (*),const SI (*), int , double, double );
void write_gadget(FILE (*), const PARTICLE (*),const SI (*));
void write_ics(FILE (*),const SI (*),int);
void write_tipsy(FILE (*), const TIPSY_STRUCTURE (*));
