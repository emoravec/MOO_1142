M2 only fit to MOO1142 using XRay starting points, with subcluster done by Jack Orlowski-Scherer in May 2026 and uploaded to Box on May 13th, 2026. 

Jack posted in Slack, saying:

My preferred model (M2 only, initialized with XRay best fit params, with subcluster) can be found here with the covariance. The covariance is now "baked in" to the model, so you can access the covariance by doing:

import dill as pk

with open("results_final_fit.dill", "rb") as f:
    metamodel, state = pk.load(f)

cov = metamodel.cov

metamodel nicely tracks everything you need to work with the covariance, for example metamodel.par_names gets you the parameter names in the same order as cov . I can also put the covariance in a text file or something like that if that's easier, but it's context free that way