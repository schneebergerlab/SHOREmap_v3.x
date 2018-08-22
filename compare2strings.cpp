#include   <string>
#include <iostream>
using namespace std;
void compare2strings(string refs, string ecos, string* diff_ecos)
{
   if(refs.size() != ecos.size()) 
   {
       cout << "strings must be the same length" << endl;
   }
   if(refs.size() == ecos.size() && ecos.size()==0)
   {
       cout << "strings are null" << endl;
   } 
   
   (*diff_ecos) = "";
   for(unsigned long istr = 0; istr < refs.size(); istr ++)
   {
       if(refs[istr] == ecos[istr]) (*diff_ecos) += "-";
       else                         (*diff_ecos) += ecos[istr];
   }
}
