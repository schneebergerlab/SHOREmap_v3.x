# Install libXm.so, libXt.so if they have not been installed  (necessary for dislin)

  --check if libXm.so* has been installed,
     
     ls /usr/lib/libXm.so*
     
  --if 'ls: cannot access /usr/lib/libXm.so*: No such file or directory' occurs, install it with
     
     sudo apt-get install libmotif4
     
  --check again (and the corresponding files should be listed),
     
     ls /usr/lib/libXm.so*
     
  --check if libXt.so*
     
     ls /usr/lib/i386-linux-gnu/libXt.so* or ls /usr/lib/x86_64-linux-gnu/libXt.so*
     
  --if (ls: cannot access /usr/lib/i386-linux-gnu/libXt.so*: No such file or directory), install it with
    
     sudo apt-get install libxt-dev
     #(answer 'y' if asked to continue!)
     
  --then check again,
     
     ls /usr/lib/i386-linux-gnu/libXt.so*
     
----------------------------------------------------------------------------------------------------
# When compiling SHOREmap

  --if '/usr/lib/ld: warning: libXm.so.3, needed by ./dislin10.3/libdislin_d.so, not found (try 
  using -rpath or -rpath-link)' occurs during compling with 'make', and you have installed 
  libmotif4, do the following:
     
     mv /usr/lib/libXm.so.4 /usr/lib/libXm.so.3
     
  --Then 'make' again, it should work.
