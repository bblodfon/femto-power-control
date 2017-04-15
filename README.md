# femtoPowerControl
My NTUA thesis code

This code is very flexible since you can randomly create FAPS (Femto Access Points) inside a Macrocell (as 
many as you want!) and lets you choose the number of macro-users (also randomly placed inside the Macro-cell). 
Also, it lets you choose the number of Non-Real-Time macro-users (data users - need high data rate services - 2.4 Mbps).
Subsequently, the rest of the macro-users are Real-Time users (voice users mostly - data rate of 64Kbps).
Each FAP has one Real-Time user and one Non-Real-Time by default.

Then, the John Nash equilibrium algorithm which is implemented, finds the equilibrium powers of all the users
(voice/data - femto and macro-users) and creates the desired graphs such as the topology of the mobile network, 
the mobile users' achieved power vs iterations of the algorithm (it is the 'k' variable in the code) needed to achieve the equilibrium, the data Rates achieved as well as the utilites of every user (both macro and femto users).

See the zompolas_femtopowercontrol pdf for a full explanation of my thesis.

The femtopowercontrol.m was tested in Matlab (version:R2011a). Run: >> femtopower
