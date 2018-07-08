# femtoPowerControl
My NTUA (National Technical University of Athens) thesis code

## Abstract

In this thesis the problem of efficient power allocation in the uplink of two-tier closed-access femtocell networks is addressed. Specifically, a single CDMA macrocell is assumed, where Ν femtocells reside within the macrocell.Within the proposed framework, which supports multiple services, appropriate utility functions are adopted to reflect users’ degree of satisfaction with respect to their actual throughput requirements and the corresponding power consumption. The overall problem is formulated as a non-cooperative game where users aim selfishly at maximizing their utility-based performance while taking into account the interference caused by both the CDMA macrocell and the neighboring femtocells.
The existence and uniqueness of a Nash equilibrium point of the proposed Multi Service Two-Tier Power Control Game with Pricing (MTTPG) is proven, at which all users have achieved a targeted SINR threshold value or transmit with their maximum power, leading essentially to an SINR-balanced system. Moreover, a distributed iterative algorithm for reaching MTTPG game’s equilibrium is provided. Finally, the operational effectiveness of the proposed approach is evaluated through modeling and simulation, while its superiority is illustrated via presenting various scenaria of the proposed framework.

### Notes
- The code is very flexible since you can randomly create FAPS (Femto Access Points) inside a Macrocell (as 
many as you want!) and lets you choose the number of macro-users (also randomly placed inside the Macro-cell). 
Also, it lets you choose the number of Non-Real-Time macro-users (data users - need high data rate services - 2.4 Mbps).
Subsequently, the rest of the macro-users are Real-Time users (voice users mostly - data rate of 64Kbps).
Each FAP has one Real-Time user and one Non-Real-Time by default.

- The John Nash equilibrium algorithm which is implemented, finds the equilibrium powers of all the users
(voice/data - femto and macro-users) and creates the desired graphs such as the topology of the mobile network, 
the mobile users' achieved power vs iterations of the algorithm (it is the 'k' variable in the code) needed to achieve the equilibrium, the data Rates achieved as well as the utilites of every user (both macro and femto users).

The femtopowercontrol.m was tested in Matlab (version: R2011a). Run: >> femtopower

The thesis report is the pdf titled: "M.Eng Thesis - Optimal Power Allocation in the uplink of Two-Tier Wireless Femtocell Networks".
