#############################################################################################
#                                                                                           #
#                                                                                           #
#############                 Fix needed nodes prior to EM                         ##########
#                                                                                           #
#                                                                                           #
#############################################################################################

# Load the HDBN
# List the nodes that should be fixed manually even before elicitation
# == some nodes are defined to reduce CPT from too many parents.
# == "Encoder node"
# Those nodes are "deterministic-ish" and should not serve EM or the network in general as "garbage solution finder"/find weird EM shortcuts
# It means "Hard coding" assumptions but its okay because they only run on logic
# e.g. If Temp is warm-ish and Chl-a is high-ish → “Above_Usual” likely
# If Temp is cool-ish and Chl-a is low-ish → “Below_Usual” likely
# Otherwise → “Usual” likely

# Here we can list easily the nodes that we want fixed, read their CPT in the shiny app and modify them manually 


# INIT ----



