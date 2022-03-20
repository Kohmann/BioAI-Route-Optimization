import kotlinx.serialization.*
import kotlinx.serialization.json.*
import java.util.Collections
import kotlin.math.max
import kotlin.math.roundToInt
import kotlin.random.Random
import kotlin.system.measureTimeMillis
import sun.misc.Signal
import sun.misc.SignalHandler


@Serializable
data class Depot(val return_time: Int, val x_coord: Int, val y_coord: Int)

@Serializable
data class Patient(val care_time: Int, val demand: Int, val end_time: Int,
                      val start_time: Int, val x_coord: Int, val y_coord: Int)

@Serializable
data class Instance(val instance_name: String, val nbr_nurses: Int,
                    val capacity_nurse: Int, val benchmark: Float,
                    val depot: Depot, val patients: Map<String, Patient>,
                    val travel_times: List<List<Float>>)

data class Data(val instance_name: String, val nbr_nurses: Int,
                    val capacity_nurse: Int, val benchmark: Float,
                    val depot: Depot, val patients: List<Patient>,
                    val travel_times: List<List<Float>>)

fun getResourceAsText(path: String): String? =
        object {}.javaClass.getResource(path)?.readText()
fun createDataclass(filename: String): Data {
    val text = getResourceAsText(filename) ?: error("Could not read file")
    val obj = Json.decodeFromString<Instance>(text)
    return Data(obj.instance_name, obj.nbr_nurses, obj.capacity_nurse,
            obj.benchmark, obj.depot, obj.patients.values.toList(), obj.travel_times)
}


class Individual(private val data: Data,
                 private val phiMax: Float = 1F,
                 initGene: List<Int>? = null) {

    private val looseMutation = true // feasible but not optimal
    private val numNurses = data.nbr_nurses
    private val numPatients = data.patients.size
    var phi: Float = if (phiMax == 1F) Random.nextFloat()
                    else if (phiMax < 1F) phiMax
                    else Random.nextInt(0, phiMax.toInt()) + Random.nextFloat()

    var gene: List<Int> = initGene ?: (IntArray(numPatients){ it } + IntArray(numNurses-1){ -(it+1) }).asList().shuffled()

    fun saveToFile() {
        val solutionGene = mutableListOf<Int>()
        gene.forEach { solutionGene.add(it+1) }
        val solution = splitListOnNegative(solutionGene, 1)
        println(solution)
    }
    fun getFitness(timePenalty: Float = 0F, route: List<Int> = gene): Float {
        return getTravelDistance(route) + timePenalty * getTimeWindowViolation(route) + getCapacityViolation(route)
    }
    private fun getTravelDistance(route: List<Int> = gene): Float {
        var travelTime = 0F
        var prevLocation = -1

        for (newLocation in route) {
            // If route is finished or not
            if (newLocation < 0) {
                // Back to depot
                travelTime += data.travel_times[prevLocation+1][0]
                prevLocation = -1 // depot
            }
            else {
                travelTime += data.travel_times[prevLocation+1][newLocation+1]
                prevLocation = newLocation
            }
        }
        travelTime += data.travel_times[prevLocation+1][0] // last step back to depot
        return travelTime
    }
    private fun getTimeWindowViolation(route: List<Int> = gene): Float {
        var prevLocation = -1
        var timeViolation = 0F
        var usedCareTime = 0F

        for (newLocation in route) {
            // If route is finished or not
            if (newLocation < 0) {
                prevLocation = -1 // depot
                usedCareTime = 0F
            }
            else {
                // add traveling time
                usedCareTime += data.travel_times[prevLocation+1][newLocation+1]

                // if arrived early,  start at start time
                if (usedCareTime < data.patients[newLocation].start_time)
                    usedCareTime = data.patients[newLocation].start_time.toFloat()

                // add caring time
                usedCareTime += data.patients[newLocation].care_time.toFloat()

                // is total time used less than end time
                timeViolation += max(0F, usedCareTime - data.patients[newLocation].end_time.toFloat())
                prevLocation = newLocation
            }
        }

        return timeViolation + max(0F, timeViolation - data.depot.return_time)
    }
    private fun getCapacityViolation(route: List<Int> = gene): Float {
        var sumDemand = 0F
        var capacitySum = 0F
        for (patient in route) {
            // If route is finished or not
            if (patient < 0) {
                sumDemand = 0F
            }
            else {
                sumDemand += data.patients[patient].demand
                if (sumDemand > data.capacity_nurse)
                    capacitySum += (sumDemand - data.capacity_nurse)
            }
        }
        return capacitySum
    }

    fun mutate(mutateProbability: Float, timePenalty: Float) {
        // chance of mutating
        if (Random.nextFloat() < mutateProbability) {
            when (Random.nextInt(0,5)) {
                0 -> swapRoutePatients(timePenalty)
                1 -> insertMutation(timePenalty)
                2 -> swapTwoBetweenRoutes(timePenalty)
                3 -> movePatientToRoute(timePenalty)
                4 -> splitRoute(timePenalty)
            }
            mutatePhi()
        }

    }
    private fun mutatePhi() {
        val randNum = smile.stat.distribution.GaussianDistribution(0.0, (phiMax * 0.1F).toDouble()).rand()
        val newPhi =  (phi.toDouble() + randNum).toFloat()

        if ( (newPhi >= 0) and (newPhi <= phiMax) ) {
                phi = newPhi
        }

    }
    private fun movePatientToRoute(timePenalty: Float) {
        val routes = splitListOnNegative(gene)
        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var bestFit = Pair(listOf(-1,-1, listOf<Int>()), 0F)

        //println("\nRoute1 $route")
        for (node in route) {
            //print("   Current Node $node")
            // finds the closest other route to given node
            var shortest = Pair(-1, Float.MAX_VALUE)
            val distances = data.travel_times[node+1]

            // Finds the closest node contained in another route
            for (pos in 1 until distances.size) {
                if (pos == node+1)
                    continue
                // distance is shorter and not in current working route
                if (distances[pos] < shortest.second && pos-1 !in route)
                    shortest = Pair(pos-1, distances[pos])
            }

            // selects the closest route
            var neighborRoute = listOf<Int>()
            for (r in routes) {
                if (shortest.first in r) {
                    neighborRoute = r
                    break
                }
            }

            if (neighborRoute.isEmpty()) {
                //println("${shortest.first} not found in $routes, \n $route")
                //error("Neighbor route is EMPTY")
                continue
            }


            val combinedFitness = getFitness(timePenalty, route) + getFitness(timePenalty, neighborRoute)
            bestFit = Pair(listOf(-1,-1, listOf<Int>()), combinedFitness)

            var tempRoute = route.toMutableList()
            tempRoute.remove(node)
            val fitnessCurrentRoute = getFitness(timePenalty, tempRoute)

            //println("   Neighbor route $neighborRoute, Fitness to beat: $combinedFitness")

            for (i in 0.. neighborRoute.size) {

                tempRoute = neighborRoute.toMutableList()
                tempRoute.add(i, node)
                //print("\t\t$tempRoute")
                val fitnessNeighbor = getFitness(timePenalty, tempRoute)
                val newCombinedFitness = fitnessNeighbor + fitnessCurrentRoute
                //print(", fitness: $newCombinedFitness")
                if (newCombinedFitness < bestFit.second) {
                    //print(",   New best fitness")
                    bestFit = Pair(listOf(node, i, neighborRoute), newCombinedFitness)
                }

            }
            if (bestFit.first[0] != -1 && looseMutation)
                break

        }
        if (bestFit.first[0] != -1) {
            val (node, index, inRoute) = bestFit.first
            if (node !is Int) error("node is not an int")
            if (index !is Int) error("index is not an int")

            val removedNodeRoute = route.toMutableList()
            removedNodeRoute.remove(node)
            routes.add(removedNodeRoute)

            var newRoute = inRoute as MutableList<Int>
            newRoute = newRoute.toMutableList()

            routes.remove(newRoute)

            newRoute.add(index, node)
            routes.add(newRoute)

        }
        else {
            routes.add(route)
        }
        gene = concatenateGene(routes)
        //validGene(gene)

    }
    private fun swapTwoBetweenRoutes(timePenalty: Float) {
        val routes = splitListOnNegative(gene)

        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var bestFit = Pair(listOf(-1,-1, listOf<Int>()), getFitness(timePenalty, route))

        //println("\nRoute1 $route")
        for (node in route) {
            //println("Current Node $node")
            // finds the closest other route to given node
            var shortest = Pair(-1, Float.MAX_VALUE)
            val distances = data.travel_times[node+1]

            // Finds the closest node contained in another route
            for (pos in 1 until distances.size) {
                if (pos == node+1)
                    continue
                // distance is shorter and not in current working route
                if (distances[pos] < shortest.second && pos-1 !in route)
                    shortest = Pair(pos-1, distances[pos])
            }

            // selects the closest route
            var neighborRoute = listOf<Int>()


            for (r in routes) {
                if (shortest.first in r) {
                    neighborRoute = r
                    break
                }
            }

            if (neighborRoute.isEmpty()) {
                //println("${shortest.first} not found in $routes, \n $route")
                //error("Neighbor route is EMPTY")
                continue
            }

            //println("Route2 $neighborRoute")

            var bestNeighborFitness = getFitness(timePenalty, neighborRoute)

            // swap and calculate fitness
            for (neighborNode in neighborRoute) {
                val temproute1 = route.toMutableList()
                val temproute2 = neighborRoute.toMutableList()

                val index1 = temproute1.indexOf(node)
                val index2 = temproute2.indexOf(neighborNode)
                temproute1.remove(node)
                temproute2.remove(neighborNode)
                temproute1.add(index1, neighborNode)
                temproute2.add(index2, node)

                val fitness1 = getFitness(timePenalty, temproute1)
                val fitness2 = getFitness(timePenalty, temproute2)
                // better for route1 and better for route2
                if (fitness1 < bestFit.second && fitness2 < bestNeighborFitness) {
                    bestNeighborFitness = fitness2
                    bestFit = Pair(listOf(node,neighborNode, neighborRoute), fitness1)
                }
            }

            //if (bestFit.first[0] != -1) {
            //    //println("Init Route1 fitness ${getFitness(timePenalty, route)}, Route2 fitness: ${getFitness(timePenalty,neighborRoute)}")
            //    //println("Best fit, Node ${bestFit.first[0]} and ${bestFit.first[1]}")
            //}

        }

        if (bestFit.first[0] != -1) {
            //println("\nBest overall: ${bestFit.first[0]} with ${bestFit.first[1]}, fitness ${bestFit.second}")

            val newRoute1 = route.toMutableList()
            val newRoute2 = bestFit.first[2] as MutableList<Int>
            routes.remove(newRoute2)
            val node1 = newRoute1.indexOf(bestFit.first[0])
            val node2 = newRoute2.indexOf(bestFit.first[1])

            newRoute1[node1] = bestFit.first[1] as Int
            newRoute2[node2] = bestFit.first[0] as Int

            routes.add(newRoute1)
            routes.add(newRoute2)

        }
        else {
            routes.add(route)
            //print("No good switch")
        }
        // route to gene and save it
        gene = concatenateGene(routes)
        //validGene(gene)
    }
    private fun swapRoutePatients(timePenalty: Float) {
        // Swap two nodes with each other in a route
        val routes = splitListOnNegative(gene)

        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        //var reducedRoute = route.toMutableList()
        var bestFit = Pair(listOf(-1,-1), getFitness(timePenalty, route))

        for (nodePos1 in route.indices) {
            for (nodePos2 in (nodePos1+1) until route.size) {
                val reducedRoute = route.toMutableList() // back to original
                Collections.swap(reducedRoute, nodePos1, nodePos2)
                val fitness = getFitness(timePenalty, reducedRoute) // TODO("Make time penalty available for Individual class")
                if (fitness < bestFit.second)
                    bestFit = Pair(listOf(nodePos1, nodePos2), fitness)

            }
            if (bestFit.first[0] != -1 && looseMutation)
                break

        }
        val newRoute = route.toMutableList()
        if (bestFit.first[0] != -1)
            Collections.swap(newRoute, bestFit.first[0], bestFit.first[1])

        routes.add(newRoute)

        // route to gene and save it
        gene = concatenateGene(routes)
        //validGene(gene)
    }
    private fun insertMutation(timePenalty: Float) {
        // positions a node in a better spot within a route
        val routes = splitListOnNegative(gene)

        // choose a random route
        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var bestFit= Pair(route, getFitness(timePenalty, route))

        for (node in route) {

            val reducedRoute = route.toMutableList()
            reducedRoute.remove(node)

            for (position in 0 .. reducedRoute.size) {
                val newRoute = reducedRoute.toMutableList() // make copy of list
                newRoute.add(position, node)
                val fitness = getFitness(timePenalty, newRoute)

                if (fitness < bestFit.second)
                    bestFit = Pair(newRoute, fitness)
            }
            if (bestFit.first[0] != -1 && looseMutation)
                break

        }
        routes.add(bestFit.first)
        // route to gene and save it
        gene = concatenateGene(routes)
        //validGene(gene)

    }
    private fun splitRoute(timePenalty: Float) {

        if (gene.last() < 0) {
            val unusedVehicle = gene.last()

            val initFitness = 0F
            var bestSplit = Pair(-1, initFitness)
            for (i in gene.indices) {
                val tempGene = gene.toMutableList()
                tempGene.remove(unusedVehicle)

                tempGene.add(i, unusedVehicle)
                val fitness = getFitness(timePenalty, tempGene)
                if (fitness < bestSplit.second)
                    bestSplit = Pair(i, fitness)
            }
            if (bestSplit.first != -1) {
                val geneCopy = gene.toMutableList()
                geneCopy.remove(gene.last())
                geneCopy.add(bestSplit.first, unusedVehicle)
                gene = geneCopy.toList()
                //println("Split Route, Old fitness: $initFitness, new Fitness: ${getFitness(10F,gene)}")
            }
        }
        //validGene(gene)
    }

    fun heuristicCrossoverSwap(father: Individual, timeSimilarity: Boolean): Individual {
        // Heuristic crossover
        // https://www.sciencedirect.com/science/article/pii/S095418100100005X?via%3Dihub

        val cut = Random.nextInt(0, data.patients.size) //

        val gene1 = gene.toMutableList()
        val gene2 = father.gene.toMutableList()

        // reorders the genes
        val tempGene1 = gene1.subList(gene1.indexOf(cut), gene1.size)
        tempGene1.addAll(gene1.subList(0, gene1.indexOf(cut)))
        val tempGene2 = gene2.subList(gene2.indexOf(cut), gene2.size)
        tempGene2.addAll(gene2.subList(0, gene2.indexOf(cut)))

        //println(tempGene1)
        //println(tempGene2)

        val newGene = mutableListOf<Int>()
        var chosenNode = cut
        newGene.add(chosenNode)

        for (i in 0 until gene.size-1) {


            if (tempGene1[i] != chosenNode) {
                Collections.swap(tempGene1, i, tempGene1.indexOf(chosenNode))
            }
            if (tempGene2[i] != chosenNode) {
                Collections.swap(tempGene2, i, tempGene2.indexOf(chosenNode))
            }

            val node1 = tempGene1[i]
            val node2 = tempGene2[i]

            val distNode1: Float
            val distNode2: Float

            if (timeSimilarity) {
                distNode1 = timeToNextPatient(tempGene1, node1)
                distNode2 = timeToNextPatient(tempGene2, node1)
            }
            else {
                distNode1 = distToNextPatient(tempGene1, node1)
                distNode2 = distToNextPatient(tempGene2, node2)
            }

            // selects the node with the smallest distance
            chosenNode = if (distNode1 < distNode2)  tempGene1[i+1] else tempGene2[i+1]
            newGene.add(chosenNode)
        }
        validGene(newGene)
        return Individual(data, father.phiMax, newGene)
    }
    private fun timeToNextPatient(gen: List<Int>, nod: Int):Float {
        val nextNode = gen[gen.indexOf(nod) + 1]
        if (nextNode < 0)
            return data.depot.return_time.toFloat()
        else
            return (data.patients[nextNode].end_time - data.patients[nextNode].care_time).toFloat()
    }
    private fun distToNextPatient(gen: List<Int>, nod: Int):Float {
        val dist: Float
        val nextNode = if (gen.indexOf(nod) != -1) gen.indexOf(nod) + 1 else error("Node not in gene")

        if (nod >= 0) {
            if (gen[nextNode] >= 0)
                dist = data.travel_times[nod+1][gen[nextNode] + 1]
            else
                dist = data.travel_times[nod+1][0]
        }
        else {
            if (gen[nextNode] >= 0)
                dist = data.travel_times[0][gen[nextNode] + 1]
            else
                dist = 0F
        }

        return dist
    }

    fun randomRoute(): List<Int> {
        // returns a random route
        val routes = splitListOnNegative(gene)
        var randRoute = routes.random()
        while (randRoute.isEmpty()) randRoute = routes.random()
        return randRoute
    }
    fun removeAndInsert(patients: List<Int>): Individual {
        // removes all given patients from gene, and reinserts them in anywhere in the gene
        // where they minimize the fitness the most
        val tempGene = gene.toMutableList()
        tempGene.removeAll(patients)

        val routes = splitListOnNegative(tempGene)

        // inserts the patients into the best possible feasible position
        // creates a new route of not possible
        for (patient in patients) {
            val insertionCosts = mutableMapOf<List<Int>, Pair<Int, Float>>()

            for (route in routes) {
                if (route.isEmpty()) continue

                // insertion cost in route
                for (i in 0..route.size) {
                    val testRoute = route.toMutableList()
                    val feasibilityError = getTravelDistance(testRoute) + getCapacityViolation(testRoute)
                    insertionCosts[route] = Pair(i, feasibilityError)
                }
            }

            if (Random.nextFloat() <= 0.8F) {
                var index = -1
                var bestRoute = mutableListOf<Int>()
                for (r in insertionCosts.keys) {
                    val idx = insertionCosts[r]?.first
                    val cost = insertionCosts[r]?.second
                    if (cost == 0F) {
                        index = idx ?: error("idx is null")
                        bestRoute = r.toMutableList()
                        break
                    }
                }
                if (index == -1) { // if no feasible
                    routes.size < data.nbr_nurses
                    routes.add(listOf(patient))
                }
                else {
                    routes.remove(bestRoute)
                    bestRoute.add(index, patient)
                    routes.add(bestRoute)
                }
            }
            else {
                if (insertionCosts.isNotEmpty()) {
                    val newRoute = insertionCosts.keys.first().toMutableList()
                    val index = insertionCosts[newRoute]?.first ?: error("Index is null")

                    routes.remove(newRoute)
                    newRoute.add(index, patient)
                    routes.add(newRoute)
                }
                else {
                    if (routes.size < data.nbr_nurses) {
                        routes.add(listOf(patient))
                    }
                    else {
                        val newRoute = routes[routes.lastIndex].toMutableList()
                        routes.remove(newRoute)
                        newRoute.add(patient)
                        routes.add(newRoute)
                    }

                }

            }

        }

        val newGene = concatenateGene(routes)
        validGene(newGene)
        return Individual(data, phiMax, newGene)
    }

    fun removeCluster(parentPhi: Float, timePenalty: Float): Individual {
        // selects a random route and find the patient with the longest arc.
        // every patient after this arc are unattached and reinserted in the best feasible way
        val routes = splitListOnNegative(gene)
        // select a random non-empty route
        var selectedRoute: List<Int>
        do {
            selectedRoute = routes.random()
        } while (selectedRoute.isEmpty())
        // remove selected route from list of all routes
        routes.remove(selectedRoute)

        var longestArc = 0F
        var longestArcIndex = 0
        var previousNode = selectedRoute.first()

        // finds the longest arc in route
        for (nextNodeIndex in 1 until selectedRoute.size ) {
            val nextNode = selectedRoute[nextNodeIndex]
            val arcDistance = data.travel_times[previousNode + 1][nextNode + 1]

            //println("$previousNode to $nextNode, distance $arcDistance")

            if (arcDistance > longestArc) {
                longestArc = arcDistance
                longestArcIndex = nextNodeIndex
            }
            previousNode = nextNode
        }

        val shortenedRoute = selectedRoute.subList(0, longestArcIndex)
        val unassignedNodes = selectedRoute.subList(longestArcIndex, selectedRoute.lastIndex+1)
        // add shortened route back in
        routes.add(shortenedRoute)

        val newGene = concatenateGene(routes).toMutableList()

        for (unassignedNode in unassignedNodes) {
            var bestFit = Pair(-1, Float.MAX_VALUE)
            for (i in 0.. newGene.size) {

                val testGene =  newGene.toMutableList()
                testGene.add(i, unassignedNode)
                val fitness = getFitness(timePenalty, testGene)
                if (fitness < bestFit.second)
                    bestFit = Pair(i, fitness)
            }
            newGene.add(bestFit.first, unassignedNode)
        }

        val newIndividual = Individual(data, phiMax ,newGene)
        newIndividual.phi = (phi + parentPhi) / 2
        validGene(newGene)
        return newIndividual

        /*
        val cluster = badRoute.slice(longestDistIdx..badRoute.lastIndex)
        badRoute.toMutableList().removeAll(cluster)

        val arcdists =  data.travel_times[cluster[0]+1]
        var shortestDist = Float.MAX_VALUE
        var patientIdx = 0
        // findest shortest arc
        for (nodeIdx in arcdists.indices) {
            val dist = arcdists[nodeIdx]
            if (dist < shortestDist && nodeIdx != cluster[0]) {
                shortestDist = dist
                patientIdx = nodeIdx
            }
        }
        var newRoute = mutableListOf<Int>()
        //println("Before loop")
        //println(patientIdx)
        //routes.forEach { println(it) }
        for (route in routes) {
            if (patientIdx in route && route.isNotEmpty()) {
                newRoute = route.toMutableList()
                patientIdx = route.indexOf(patientIdx)
                //println("Updated patientID")
                break
            }

        }

        tempRoutes.remove(newRoute.toList())

        if (newRoute.isNotEmpty()) {
            newRoute.addAll(patientIdx+1, cluster)
            tempRoutes.add(newRoute.toList())

            val newGene = mutableListOf<Int>()

            var nurse = -1
            for (route in routes) {
                if (route.isNotEmpty()) {
                    route.forEach { newGene.add(it)}
                    newGene.add(nurse--)
                }
            }
            gene = newGene
        }
        else
            println("Error")

         */
    }

    fun insertionHeuristicSimple(used: List<Int>) {
        // Greedy construction heristic that selects the closest and has the earliest
        // start time. If demand exceeds capacity, a new route is created.
        // patientStartTime = [ (index, patient), ... ], index = patientID aka node number
        val patientList = data.patients.mapIndexed { index: Int, patient -> index to patient  }.toMutableList()

        var usedNurses = 1
        val initRoutes = mutableListOf<List<Int>>()
        var firstRun = true

        while (patientList.isNotEmpty() && (usedNurses < (data.nbr_nurses-2) )) {

            val patientListSorted = patientList.sortedBy { it.second.start_time + it.second.care_time }.toMutableList()
            var nodePair = patientListSorted.first()

            // assures that the first route is different
            if (firstRun && used.isNotEmpty()) {
                var i = 1
                while (nodePair.first in used) {
                    nodePair = patientListSorted[i++]
                }
                firstRun = false
            }
            patientList.remove(nodePair)

            // init new route
            val route = mutableListOf<Int>(nodePair.first)

            while (patientList.isNotEmpty()) {
                val patientMap = patientList.toMap()
                var lowestFitness = Pair(-1, Float.MAX_VALUE) // (index, fitness)
                var nrFeasibleNodes = 0


                //println("\nRoute $route, Distance from node ${route.last()}")
                val shortestDists = mutableListOf<Pair<Int,Float>>()

                for (i in 0 until data.travel_times[0].size-1) {
                    if (route.last()+1 != i && i in patientMap.keys) {
                        shortestDists.add(Pair(i, data.travel_times[route.last()+1][i]) )
                    }
                }
                shortestDists.sortBy { it.second }


                val shortestTime = mutableListOf<Pair<Int,Int>>()
                for (x in shortestDists) {
                    patientMap[x.first]?.let { Pair(x.first, it.start_time) }?.let { shortestTime.add(it) }
                }
                shortestTime.sortBy { it.second }


                for (node in shortestDists) {
                    // tests every position in route

                    val temproute = route.toMutableList()
                    temproute.add(node.first)

                    //println("$temproute, timePenalty: ${getTimeWindowViolation(temproute)}, Capacaty: ${getCapacityViolation(temproute)}")

                    // Node is feasible
                    if (getCapacityViolation(temproute) == 0F && getTimeWindowViolation(temproute) == 0F) {

                        val fitness = 10 * getTravelDistance(temproute) + (patientMap[node.first]?.start_time ?: 0)

                        if (fitness < lowestFitness.second) {
                            lowestFitness = Pair(node.first, fitness)
                        }

                        //lowestFitness = Pair(node.first, 0F)
                        nrFeasibleNodes++
                        //break
                    }

                }
                // if no feasible nodes exist
                if (nrFeasibleNodes == 0) {
                    break
                }

                // adds the best node to route
                route.add(lowestFitness.first)
                patientList.removeIf { it.first == lowestFitness.first }

            }
            usedNurses++
            initRoutes.add(route)

        }

        val route = mutableListOf<Int>()
        patientList.forEach { route.add(it.first) }
        initRoutes.add(route)
        gene = concatenateGene(initRoutes)
        //validGene(gene)
    }
    fun insertionHeuristic(used: List<Int>, backwards: Boolean = false, randChoise: Boolean = true) {
        // Very greedy construction heuristic that selects the closest node from Depot and always
        //  adds best patient with respect to distance and time window
        // patientStartTime = [ (index, patient), ... ], index = patientID aka node number
        val patientList = data.patients.mapIndexed { index: Int, patient -> index to patient  }.toMutableList()

        var usedNurses = 1
        val initRoutes = mutableListOf<List<Int>>()
        var firstRun = true


        while (patientList.isNotEmpty() && (usedNurses < (data.nbr_nurses-2) )) {

            val patientListSorted = if (backwards)
                    patientList.sortedBy { it.second.end_time - it.second.care_time }.toMutableList()
                    else
                    patientList.sortedBy { it.second.start_time + it.second.care_time }.toMutableList()

            // reverse the list or not, is random
            if (Random.nextFloat() > 0.5F)
                patientListSorted.reverse()

            var nodePair = if (randChoise) patientListSorted.random() else patientListSorted.first()

            // assures that the first route is different
            if (firstRun && used.isNotEmpty()) {
                var i = 1
                while (nodePair.first in used) {
                    nodePair = patientListSorted[i++]
                }
                firstRun = false
            }


            patientList.remove(nodePair)
            // init new route
            val route = mutableListOf<Int>(nodePair.first)
            //println("first route: $route Start time: ${nodePair.second.start_time}")

            while (patientList.isNotEmpty()) {

                var lowestFitness = Pair(0, Float.MAX_VALUE) // (index, fitness)
                var bestPatient = patientList[0].second // random init
                var nodeposition = 0
                var foundNonode = false // if a possible node exists
                var nrFeasibleNodes = 0

                for (node in patientList) {
                    // tests every position in route
                    for (position in 0..route.size) {

                        val temproute = route.toMutableList()
                        temproute.add(position, node.first)
                        //println("$temproute, travelDist: ${getTravelDistance(temproute)}, timePenalty: ${getTimeWindowViolation(temproute)}, Capacaty: ${getCapacityViolation(temproute)}")

                        if (getCapacityViolation(temproute) == 0F && getTimeWindowViolation(temproute) == 0F) { //

                            foundNonode = true
                            val fitness = if (backwards)
                                    10 * getTravelDistance(temproute) - (node.second.end_time - node.second.care_time)
                                    else
                                    10 * getTravelDistance(temproute) + node.second.start_time
                            //println("travelDist: ${getTravelDistance(temproute)}, timePenalty: ${getTimeWindowViolation(temproute)}, Capacaty: ${getCapacityViolation(temproute)} , route $temproute")
                            if (fitness < lowestFitness.second) {
                                lowestFitness = Pair(node.first, fitness)
                                bestPatient = node.second
                                nodeposition = position
                            }
                            nrFeasibleNodes++
                        }
                    }
                }
                if (!foundNonode) {
                    break
                } // if no feasible nodes exist

                // adds the best node to route
                route.add(nodeposition, lowestFitness.first)
                patientList.remove(Pair(lowestFitness.first, bestPatient))
            }
            usedNurses++
            //if (backwards)
            //    route.reverse()
            initRoutes.add(route)
        }
        // add remaining un-assigned patients
        val route = mutableListOf<Int>()
        patientList.forEach { route.add(it.first) }
        //if (backwards)
        //    route.reverse()
        initRoutes.add(route)

        gene = concatenateGene(initRoutes) // save gene
        //validGene(gene)
    }
    fun insertionHeuristicMultible() {
        // Selects random patients and adds on the nearest patient in a round-robin fashion
        val numberOfRoutes = Random.nextInt(1, data.nbr_nurses+1)

        val routes = mutableListOf<MutableList<Int>>()
        val availablePatients = MutableList(data.patients.size){it}.shuffled().toMutableList()
        for (i in 0 until numberOfRoutes) {
            routes.add(mutableListOf(availablePatients.removeAt(0)))
        }

        while (availablePatients.isNotEmpty()) {
            for (route in routes){
                if (availablePatients.isEmpty()) break
                val prevNode = route.last()
                // closest node from prevNode
                var shortestDists = Pair(-1,Float.MAX_VALUE)
                for (i in availablePatients) {
                    val dist = data.travel_times[prevNode+1][i+1]
                    if (dist < shortestDists.second)
                        shortestDists = Pair(i, dist)
                }
                route.add(shortestDists.first)
                availablePatients.remove(shortestDists.first)
            }
        }
        gene = concatenateGene(routes)
        //validGene(gene)
    }

    private fun splitListOnNegative(arr: List<Int>, onValue: Int = 0): MutableList<List<Int>> {
        // Creates a list of lists (graph) of chromosome representation and returns the graph
        var startIndex = 0
        val splitedList: MutableList<List<Int>> = ArrayList()
        val arrLength = arr.size
        for (i in 0 until arrLength) {
            if (arr[i] < onValue) {
                splitedList.add(arr.subList(startIndex, i))
                startIndex = i+1
            }
        }
        splitedList.add(arr.subList(startIndex, arrLength))
        return splitedList
    }
    private fun validGene(arr: List<Int>) {
        for (i in 0..99) {
            if (i !in arr)
                error("Not a valid gene, Patient problem: ${arr.size}, $arr")
        }
        for (i in -1 downTo -24) {
            if (i !in arr)
                error("Not a valid gene. Nurse problem: ${arr.size}, $arr")
        }
        if (100 in arr)
            error("Not a valid gene: includes 100, $arr")
        if (arr.size != 124)
            error("Not a valid gene, size error: ${arr.size}, $arr")

    }
    private fun concatenateGene(routes: List<List<Int>>): List<Int> {
        // concatenates a sequence of graphs to chromosome with negative values as delimiters
        val newGene = mutableListOf<Int>()

        var nurse = -1
        for (route in routes) {
            if (route.isNotEmpty()) {
                //route.forEach { newGene.add(it)}
                newGene.addAll(route)
            }
            if (-nurse < data.nbr_nurses)
                newGene.add(nurse--)
        }

        //adds remaining nurses to gene
        val nursesLeft = data.nbr_nurses + nurse
        repeat(nursesLeft) {
            if (nurse == -data.nbr_nurses) print("GMMM")
            newGene.add(nurse--)
        }

        return newGene
    }

    fun writeSolution() {
        // Writes solution to console
        val routes = splitListOnNegative(gene)
        println("\nNurse capacity: ${data.capacity_nurse}")
        println("Depot return time: ${data.depot.return_time}")
        var nurseNr = 1
        for (route in routes) {
            if (route.isEmpty()) continue
            print("Nurse %-2d   %-7.2f".format(nurseNr++, getTravelDistance(route)))
            print("  D(0) -> ")


            var prevPatient = -1
            var endTime = 0F
            for (patient in route) {
                val travelDist = data.travel_times[prevPatient+1][patient+1]
                val startTime = data.patients[patient].start_time
                val careTime = data.patients[patient].care_time
                endTime += travelDist
                if (endTime < startTime)
                    endTime = startTime.toFloat()


                print("%-3d (%-6.2f-%-6.2f)".format(patient+1, endTime, endTime + careTime) +
                        "[%-4d-%-4d]".format(data.patients[patient].start_time,
                                    data.patients[patient].end_time))
                endTime += careTime
                print(" \t -> ")
                prevPatient = patient
            }
            print("D($endTime) \n")
        }
        println("Objective value (total duration): %.2f".format(getTravelDistance(gene)))
    }
}



class GeneticAlgorithm(private val data: Data,
                       private val param: GAparameters) {

    private val sizeOfPopulation = param.sizeOfPopulation
    private val tournamentSize = param.tournamentSize
    private val mutateProbability = param.mutateProbability
    private var timePenalty = param.timePenalty
    private val randomSelection = param.randomSelection
    private val crowdingSelectionProb = param.crowdingSelectionProb
    private val crossoverProb = param.crossoverProb
    private val diversityThreshold = param.diversityThreshold
    private val useConstructionHeuristic = param.constructionHeurstics


    var population = mutableListOf<Individual>()
    lateinit var fittestIndividual: Individual

    private fun uniquePhenoTypes(): Float {
        val fitness = mutableListOf<Float>()
        population.forEach { fitness.add(it.getFitness()) }

        return fitness.distinct().size / fitness.size.toFloat()
    }

    private fun initPopulation() {
        if (useConstructionHeuristic)
            population.addAll(constructionHeuristic(sizeOfPopulation))
        else
            repeat(sizeOfPopulation) { population.add(Individual(this.data, param.phiMax)) }
    }
    private fun constructionHeuristic(num: Int): List<Individual> {
        val pop = mutableListOf<Individual>()
        repeat(num) {
            val indiv = Individual(this.data, param.phiMax)
            indiv.insertionHeuristicMultible()
            pop.add(indiv)
        }
        return pop

        /*
        var used = mutableListOf<Int>()

        var third = (sizeOfPopulation / 3)
        repeat(third) {
            val indiv = Individual(this.data, param.phiMax)
            indiv.insertionHeuristic(used, listOf<Boolean>(true, false).random(), false)
            used.add(indiv.gene[0])
            population.add(indiv)
        }

        var adjust = if (third >= 100) 99 else third

        used = mutableListOf<Int>()
        // if half is odd
        repeat(adjust) {
            val indiv = Individual(this.data, param.phiMax)
            indiv.insertionHeuristicSimple(used)
            used.add(indiv.gene[0])
            population.add(indiv)
        }

        repeat(sizeOfPopulation - third - adjust) {
            population.add(Individual(this.data, param.phiMax))
        }
         */
    }

    private fun getAverageFitness(): Double {
        val fitness = mutableListOf<Float>()
        population.forEach { fitness.add(it.getFitness(timePenalty)) }
        val fittestIndex = fitness.indices.minByOrNull{ fitness[it] } ?: error("This should not run")
        fittestIndividual = population[fittestIndex]
        //dataStorage.fittest.add(population[fittestIndex])
        return fitness.average()
    }
    private fun getNBestIndividuals(num: Int): List<Individual> {
        val fitness = mutableListOf<Pair<Individual,Float>>()
        population.forEach { fitness.add( Pair(it, it.getFitness(timePenalty)) ) }
        val fittestPairs = fitness.sortedBy { it.second }
        val fittest = mutableListOf<Individual>()
        for (i in 0 until num)
            fittest.add(fittestPairs[i].first)
        return fittest
    }

    fun fitMultible(numGenerations: Int, threads: Int = 1) {
        val bestIndividuals = mutableListOf<Individual>()

        repeat(threads) {
            println("NEW RUN")
            population.clear()
            this.initPopulation()
            fit(numGenerations)
            bestIndividuals.addAll(getNBestIndividuals((sizeOfPopulation / threads)))
        }
        if (threads > 1) {
            population.clear()
            population = bestIndividuals
            println("\nFINAL RUN")
            fit(numGenerations/3)
        }

    }
    private fun fit(numGenerations: Int)  {
        //this.initPopulation()
        println("Population of ${population.size} initiated")

        val printNr = if (numGenerations < 30) 1 else (numGenerations * 0.01 + 1).roundToInt()
        for (i in 1..numGenerations) {
            doGeneration()
            if (i % printNr == 0) {
                //storeCurrentState()
                print("Gen: %-3d, avgFitness: %-10.2f   Fittest: %.2f, Diversity: %.3f, phi: %.2f"
                        .format(i, getAverageFitness(), fittestIndividual.getFitness(timePenalty), uniquePhenoTypes(), fittestIndividual.phi))
                val feasible = fittestIndividual.getFitness() != fittestIndividual.getFitness(timePenalty)
                if (feasible) println(", not Feasible") else println()

            }

            if ( i == (numGenerations * 0.8).toInt() ) {
                println("Super mutating!!")
                val newpop = superMutation(population).toMutableList()
                population.clear()
                population = newpop
            }
        }
        val finalFittest = superMutation(listOf(fittestIndividual)).toMutableList().first()
        if (finalFittest.getFitness(100F) < fittestIndividual.getFitness(100F))
            fittestIndividual = finalFittest
    }
    private fun doGeneration() {
        val offsprings = mutableListOf<Individual>()

        do {
            val (mother, father) = parentSelection()
            val (offspring1, offspring2) = getOffspring(mother, father) // lag sexy avkom
            offsprings.add(offspring1)
            offsprings.add(offspring2)

        } while (offsprings.size < population.size)


        val newPopulation = survivorSelection(offsprings.toList())

        population.clear()
        population = newPopulation.toMutableList()

        if (uniquePhenoTypes() < diversityThreshold) {
            //println("Added more diversity")
            population = nFittest(newPopulation).toMutableList()
            repeat(30) {
                population.forEach { it.mutate(1F, timePenalty) }
            }
        }

    }

    private fun parentSelection(): Array<Individual> {
        val mother: Individual = getParent()
        var father: Individual = getParent()
        while (mother == father) {  father = getParent() }
        return arrayOf(mother, father)
    }
    private fun getParent(): Individual {
        // TODO("Add more advanced parent selection methods")
        if (randomSelection)
            return population.random()
        else
            return tournamentSelection()
    }
    private fun tournamentSelection(): Individual {

        var bestParent = population.random()
        var competingParent = population.random()
        repeat(tournamentSize) {
            while (bestParent == competingParent) { competingParent = population.random() }
            if (bestParent.getFitness(timePenalty) > competingParent.getFitness(timePenalty)) { // minimizing fitness
                bestParent = competingParent
            }
        }
        return bestParent
    }
    private fun getOffspring(mother: Individual, father: Individual): List<Individual> {
        val (firstChild, secondChild) = doCrossover(mother, father)
        firstChild.mutate(mutateProbability, timePenalty)
        secondChild.mutate(mutateProbability, timePenalty)
        return listOf(firstChild, secondChild)

    }

    private fun doCrossover(mother: Individual, father: Individual): List<Individual> {
        // TODO("Add more advanced crossover methods")

        if (Random.nextFloat() > crossoverProb ) {
            val firstChild = mother.removeCluster(father.phi, timePenalty)
            val secondChild = father.removeCluster(mother.phi, timePenalty)
            return listOf(firstChild, secondChild)
        }
        else {
            if (Random.nextFloat() < 0.5F)
                return heuristicMergeCrossover(mother, father)
            else
                return bestCostRouteCrossover(mother, father)
        }

        //val firstChild = alternatingEdgesCrossover(mother, father)
        //val secondChild = alternatingEdgesCrossover(father, mother)
        //return listOf(firstChild, secondChild)

    }
    private fun heuristicMergeCrossover(mother: Individual, father: Individual): List<Individual> {
        val offspring1 = mother.heuristicCrossoverSwap(father, false)
        val offspring2 = father.heuristicCrossoverSwap(mother, true)
        return listOf<Individual>(offspring1, offspring2)
    }
    private fun bestCostRouteCrossover(mother: Individual, father: Individual): List<Individual> {

        val patientsMother = mother.randomRoute()
        val patientsFather = father.randomRoute()
        val offspring1 = mother.removeAndInsert(patientsFather)
        val offspring2 = father.removeAndInsert(patientsMother)

        return listOf(offspring1, offspring2)
    }
    private fun alternatingEdgesCrossover(mother: Individual, father: Individual): Individual {
        val motherGene = mother.gene//toGraph(mother.gene)
        val fatherGene = father.gene//toGraph(father.gene)


        val childGene = mutableListOf<Int>(motherGene[0], motherGene[1])
        val unusedEdges = motherGene.toMutableList()
        unusedEdges.removeAll{ it in childGene }
        var arcGene: Int
        var arcIdx: Int
        var nextArc: Int

        for (i in 1 until motherGene.size-1) {
            arcGene = childGene[i]
            if (i % 2 == 0) {
                arcIdx = motherGene.indexOf(arcGene)
                arcIdx = if (arcIdx+1 == motherGene.size) 0 else arcIdx +1
                nextArc = motherGene[arcIdx]
            }
            else {
                arcIdx = fatherGene.indexOf(arcGene)
                arcIdx = if (arcIdx+1 == fatherGene.size) 0 else arcIdx +1
                nextArc = fatherGene[arcIdx]
            }
            // already in new gene
            if (nextArc in childGene) {
                nextArc = unusedEdges.random()
            }
            childGene.add(nextArc)
            unusedEdges.remove(nextArc)
            if (unusedEdges.isEmpty()) break
        }
        if (unusedEdges.isNotEmpty())
            childGene.addAll(unusedEdges)

        //val newGene = simpleSplitToRoutes(childGene)
        //val child = Individual(data, childGene)
        //child.addDelimiters() // adds de-limiters (nurses)


        return Individual(data, (mother.phi + father.phi)/2, childGene)
    }

    private fun survivorSelection(offsprings: List<Individual>): List<Individual> {
        // TODO("Add adjustable parameter for survivor selection")
        if (Random.nextFloat() < crowdingSelectionProb)
            return generalizedCrowding(offsprings)
        else
            return nFittest(offsprings)
    }
    private fun nFittest(offsprings: List<Individual>): List<Individual> {

        val competingPopulation = offsprings + population

        val fitnessList = competingPopulation.associateBy({it.getFitness(timePenalty)}, {it})

        val fittestIndividualsList = fitnessList.toSortedMap().values.toList()
        val numberOfIndividuals = (sizeOfPopulation * 0.4F).toInt()
        val fittestIndividuals = fittestIndividualsList.subList(0, minOf(numberOfIndividuals, fittestIndividualsList.size)).toMutableList()

        if (useConstructionHeuristic)
            fittestIndividuals.addAll(constructionHeuristic(sizeOfPopulation - fittestIndividuals.size))
        else {
            repeat(sizeOfPopulation - fittestIndividuals.size) {
                fittestIndividuals.add(Individual(this.data, param.phiMax))
            }
        }

        if (Random.nextFloat() < 0.005F) {
            println("Super mutating!!")
            return superMutation(fittestIndividuals)
        }
        return fittestIndividuals
    }
    private fun superMutation(pop: List<Individual>): List<Individual> {
        repeat(500) {
            pop.forEach { it.mutate(1F, timePenalty) }
        }
        return pop
    }
    private fun generalizedCrowding(offsprings: List<Individual>): List<Individual> {
        val newGeneration = mutableListOf<Individual>()
        for ((child, parent) in offsprings zip population) {
            val F_c = child.getFitness(timePenalty)
            val F_p = parent.getFitness(timePenalty)

            val p_c_bigger = F_c / (F_c + parent.phi * F_p) // probability for child selection if child fitness is better
            val p_c_smaller = child.phi * F_c / (child.phi * F_c + F_p) // probability for child selection if parent fitness is better
            val choice: Individual


            if (F_c > F_p) { // child has higher fitness than parent
                choice = if (Random.nextFloat() <= p_c_bigger)
                                parent
                            else
                                child
            }
            else if (F_c == F_p)
                choice = if (Random.nextFloat() <= 0.5F)
                            parent
                        else
                            child
            else {
                choice = if (Random.nextFloat() <= p_c_smaller)
                    parent
                else
                    child

            }
            newGeneration.add(choice)
        }

        return newGeneration
    }


}

data class GAparameters(
        val sizeOfPopulation: Int = 200, // 100,200, 400 works well
        val tournamentSize: Int = 4, // 4 works well
        val mutateProbability: Float = 0.3F, // 0.8G is good
        val phiMax: Float = 0.8F,
        val timePenalty: Float = 40F, // 40F is good
        val randomSelection: Boolean = false,
        val constructionHeurstics: Boolean = false, // false
        val crowdingSelectionProb: Float = 1F, // crowding og eliteism
        val crossoverProb: Float = 0.4F, // 0.1F and 0.6F is good,  removeCluster or crossover heuristics
        val diversityThreshold: Float = 0.7F // 0.8F
)


fun main() {

    // Ensures that every run is different
    smile.math.MathEx.setSeed()
    val filename = "train_6.json"
    val data = createDataclass(filename)
    val model = GeneticAlgorithm(data, GAparameters())

    // prints out the best solution when stopping GA mid-run
    Signal.handle(Signal("INT"), object : SignalHandler {
        override fun handle(sig: Signal) {
            println("\nThe fittest individual: ${model.fittestIndividual.getFitness(0F)}")
            println("\nThe fittest individual with constraint: ${model.fittestIndividual.getFitness(10F)}")
            model.fittestIndividual.saveToFile()
            //model.fittestIndividual.writeSolution()
            System.exit(0)
        }
    })


    val testing = false
    if (!testing) {
        val timeInMillis = measureTimeMillis {
            //model.fit(500)
            model.fitMultible(180, 3)
        }

        println("(The operation took $timeInMillis ms)")
        println("\nThe fittest individual: ${model.fittestIndividual.getFitness(0F)}")
        println("\nThe fittest individual with constraint: ${model.fittestIndividual.getFitness(10F)}")
        model.fittestIndividual.saveToFile()

        //model.fittestIndividual.writeSolution()
    }




}
