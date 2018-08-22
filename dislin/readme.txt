1. FLOAT-precision dislin.h is not precise enough in our projects. So DOUBLE-precision dislin_d.h 
   is used in SHOREmap. However, from the folder where dislin is installed, we cannot find this 
   header file.

   To find and use dislin_d.h, we can uncompress the installation tar file: 
   dislin-11.0.linux.i586_64.tar.gz. It is in the subfolder 'examples'. 
   Copy it to the intalled folder of dislin, where dislin.h exists.

   Now we can compile SHOREmap.

2. Here, we directly copy dislin_d.h under a subfolder of SHOREmap, together with the shared libary 
   libdislin_h.so*, which can be copied from the folder where dislin is installed. 
   
   *If dislin is installed with 'sudo dpkg -i dislin-11.0.linux.i586.deb', there is no further steps
    required for setting dislin.
   
   *If disin is installed from 'dislin-11.0.linux.i386.tar.gz', further steps are necessary:
   
   --After installation, copy these files to folder path/to/SHOREmap_v3.0/dislin
         
         cp installed/path/to/dislin/*dislin_d.* path/to/SHOREmap_v3.0/dislin
         (#answer 'y' if asked to overwrite!)
     
   --Open /etc/profile 
     
         vi /etc/profile
     
   --Press key 'esc' then key 'i' to enter the edit mode of vi,
   ----if LD_LIBRARY_PATH already exists,insert the following statement into your /etc/profile, 
         
         export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/SHOREmap_v3.0/dislin
   
   ----otherwise,
         
         export LD_LIBRARY_PATH=/path/to/SHOREmap_v3.0/dislin
         
   ----quit and save modification in /etc/profile with :wq, and source it with
         
         source /etc/profile
         
## for static compiling: 

cp /projects/dep_coupland/grp_nordstrom/privatelibs/dislin/lib/dislin_d-11.0.a SHOREmap_v3.0/dislin     
