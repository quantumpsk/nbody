# nbody
An n-body system evolution simulation

What follows here is a day by day log of the growth of this bit of code.

Saturday, April Fools Day, 01/04
I'm gonna try something here with the n body problem..
I'll start with simulating a simple solar system..
or maybe an even more simple three body system..
we shall later try to scale it for any physics problem..
lets begin with the basic class structures..

Sunday, 02/04
So ive got the basic structure down, going ahead with calculating some
of the variables..
realized i have to store all vectors as ndarrays using numpy since i cant
assign lists to a scope of a ndarray.. so thats next.

Tuesday, 04/04
Well the System and its Objects are set up and i have the complete initial state.
Next is the evolution of the system..with a time step..lets brute this shit..
well that was simple, a time.step, a calculate.state in a system.evolve and
we have what we're after.
Lets get started with plotting it..(after tonights date)

Wednesday, 05/04
So after date night(thank you, alisha!), I've got to pin down the exact way im storing all the
positions of the different objects over the course of the sim to be able to plot
the full path.. and apparently the best way to dynamically grow a numpy array
is to not do it, rather grow a list and then create the array from the list ! (courtesy stackoverflow)

Thursday, 06/04
Ok.. its been over 24 hours now that im stuck with what seems should be an easy
problem to fix..
The issue is i need to store all the positions of the different objects. I thought
of growing a list as each new position is calculated. But for some reason each
log of each object is storing not only its own position but also that of other
objects and i dont know why that is..weird...i think my compiler is possessed.
...
So for now im temporarily giving up on this method of arrays bull to store the
positions, instead im storing it in a file, maybe itll be easier.

Friday, 07/04
early in the morning, here i am sitting for one more attempt at a position
log array, to convert to np.array, to plot. I dont know yet why append is adding
more than it should, but it is, so lets work with and maybe i can learn something..
...
oh man, marathon session of coding. i solved the log and plot problem. i still
dont know why append is behaving the way it is, but some ghetto work around and
i have the positions and now the plot.
Only then did i realize i hadnt accounted for the zero-distance error where
after being attracted close to each other objects then experience immense g forces
which end up scattering them.
(This was when the ion version was forked away)

Sunday, 09/04
Yesterday was a break day, and rightfully so given I was able to obtain my first
simulation plot of a 5 body system for 40k steps. while there is still the un-
explained working of the position log list append mechanism, its been ghetto-fixed
and is workable.
Im not completely happy with the plot in itself, as in the way the physics of
the interaction of the objs plays out, i think theres something im missing in there.
For now, though, the next task is going to be fine tuning the plotting, and see
about implementing a cProfile to see where i can optimize the process the most.

Monday, 24/04
Yeah its been a while and almost as per usual i never bothered to return to this
little pet project. but thanks to my dad, here i am. I never got around to doing
the fine tuning or the cProfile bits. not gonna pursue that now, want to right
the issues with the physics of the model. something basic im getting wrong. 
lets see...
