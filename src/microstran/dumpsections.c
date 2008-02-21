/**@FILE
	little program to demonstrate that we can meaninfully parse the 'properties.txt' file
*/

#include "sections.h"
#include "sectionsparser.h"

#include <stdlib.h>

int main(int argc, char **argv){
	const char *filename = "src/microstran/properties.txt";
	if(argc==2){
		filename = argv[1];
	}

	FILE *f;
	f = fopen(filename,"r");
	if(f==NULL){
		fprintf(stderr,"Unable to open section library '%s'",filename);
		exit(1);
	}

	parse *p;
	p = parseCreateFileName(filename);

	section_library *l = NULL;
	l = section_library_create();
	parseSections(p,l);

	const char *sn = "100X100X9.0SHS";
	fprintf(stderr,"\n\nLooking up section '%s'...\n",sn);
	const section *s;
	s = section_find(l,sn);
	if(s==NULL){
		fprintf(stderr,"SECTION NOT FOUND\n");
		return 1;
	}else{
		fprintf(stderr,"SECTION FOUND\n");
		if(section_is_chs(s)){
			fprintf(stderr,"CHS\n");
			fprintf(stderr,"outside diameter = %f\n", section_chs_outside_diameter(s));
			fprintf(stderr,"thickness = %f\n\n", section_chs_thickness(s));
		}
	}
	return 0;
}

