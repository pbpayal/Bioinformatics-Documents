



sudo adduser Payal
sudo passwd ****


$ sudo rstudio-server stop
$ sudo rstudio-server start
$ sudo rstudio-server restart

To list all currently active sessions:
$ sudo rstudio-server active-sessions

To suspend an individual session:
$ sudo rstudio-server suspend-session <pid>

To suspend all running sessions:
$ sudo rstudio-server suspend-all

The suspend commands also have a "force" variation which will send an interrupt to to the session to request the termination of any running R command:
$ sudo rstudio-server force-suspend-session <pid>
$ sudo rstudio-server force-suspend-all

Taking the Server Offline

If you need to perform system maintenance and want users to receive a friendly message indicating the server is offline you can issue the following command:
$ sudo rstudio-server offline
When the server is once again available you should issue this command:
$ sudo rstudio-server online



