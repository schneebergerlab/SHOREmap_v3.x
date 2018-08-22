////// test if a const char* is composed of only '0' to '9', '.' //////
bool is_number(const char* opt)
{
   int i1 = 0;
   while (*(opt+i1))
   {
       if(!(*(opt+i1) >= '0' && *(opt+i1) <= '9' || *(opt+i1) == '.')) 
           return false;
       i1 ++;
   }
   return true;
}
