/* input:  file contains chrmosome ids and postions of reference errors             */
/* output: <"chr_id.#.ref_error_position", true/false> in global variable REFERROR  */
/* date:   2013-03-05                                                               */
/* file format: chr_id	ref_error_position
		1	1590
		2	1592
		3	1593
		4	1594
*/
bool read_referror(char* freferror);
