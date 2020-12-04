# Load distance matrix

X <- read.table("X.csv")

# Compute persistence

persistence <- TDA::ripsDiag(X=X, maxdimension=1, maxscale=Inf, dist="arbitrary", location=TRUE, library="Dionysus")

# Save results

write.table(persistence$diagram, "diagram.csv")
write(jsonlite::toJSON(persistence$cycleLocation), "cycleLocation.json")
