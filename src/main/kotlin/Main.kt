import kotlinx.serialization.*
import kotlinx.serialization.json.*
import smile.clustering.*
import java.util.Collections
import kotlin.math.max
import kotlin.math.roundToInt
import kotlin.random.Random
import kotlin.system.measureTimeMillis


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


class LinearBoundedSplit(private val data: Data, var graphRep: MutableList<Int>) {

    var graph = graphRep  // graph representation

    var pred: MutableList<MutableList<Int>> = mutableListOf()
    var potential: MutableList<MutableList<Float>> = mutableListOf()

    val sumDistance: MutableList<Float> = mutableListOf()
    val sumLoad: MutableList<Int> = mutableListOf()


    private fun initData() {
        for (nurse in 0..data.nbr_nurses ) {
            val tempDbl = MutableList(graph.size+1){10000000F}
            val tempInt = MutableList(graph.size+1){0}
            potential.add(tempDbl)
            pred.add(tempInt)
        }
        potential[0][0] = 0F
        distanceData()
    }
    private fun distanceData() {
        var travelTime = 0F
        sumDistance.add(travelTime)
        for (node in 0 until  (graph.size-1) ) {
            travelTime += data.travel_times[graph[node]][graph[node+1]]
            sumDistance.add(travelTime)
        }
        travelTime += data.travel_times[0][graph.last()]
        sumDistance.add(travelTime)

        sumLoad.add(0)
        for (node in graph) {
            sumLoad.add(sumLoad.last() + data.patients[node].demand)
        }
    }

    fun findDelimiters(): List<Int> {
        initData()
        val queue: MutableList<Int> = mutableListOf()

        for (k in 0 until data.nbr_nurses) {
            queue.clear()
            queue.add(k)

            for (i in k+1.. graph.size) {
                potential[k+1][i] = propagate(queue.first(), i, k)
                pred[k+1][i] = queue.first() //graph[queue.first()] //

                if (i < graph.size) {
                    if (!dominates(queue.last(), i, k)) {
                        while (queue.isNotEmpty() && dominatesRight(queue.last(), i, k)) {
                            queue.removeLast()
                        }
                        queue.add(i)
                    }
                    while (queue.isNotEmpty() &&
                            (sumLoad[i+1] - sumLoad[queue.first()]) > (data.capacity_nurse + 0.0001F) ) {
                        queue.removeFirst()
                    }
                }

                if (queue.isEmpty()) { break }
            }
        }

        return splitGraphToGene()
    }

    private fun splitGraphToGene(): List<Int> {
        var minCost = 1000000F
        var minIndex: Int = -1
        for (i in potential.indices) {
            if (potential[i].last() < minCost) {
                minCost = potential[i].last()
                minIndex = i
            }
            //println("Cost:${split.potential[i].last() }, minCost:$minCost")
        }
        val newGene = graph
        val newGeneTest = MutableList<Int>(minIndex){0}
        var cour = data.patients.size
        for (i in minIndex-1 downTo 0) {
            //println("i+1: ${i+1}, cour: $cour" + " value: ${pred[i+1][cour]}")
            cour = pred[i+1][cour]
            newGeneTest[i] = cour+1
            newGene.add(cour+1, -1)
        }

        println("Best k: $minIndex, List:$newGeneTest")
        return newGene
    }
    private fun propagate(i: Int, j: Int, k: Int): Float {
        return potential[k][i] + sumDistance[j] - sumDistance[i+1] +
                data.travel_times[0][graph[i]] + data.travel_times[0][graph[j-1]]
    }
    private fun dominates(i: Int, j: Int, k: Int): Boolean {

        return sumLoad[i] == sumLoad[j] &&
                ((potential[k][j] + data.travel_times[0][graph[j]]) >
                (potential[k][i] + data.travel_times[0][graph[i]] +
                        sumDistance[j+1] - sumDistance[i+1] - 0.0001F))
    }
    private fun dominatesRight(i: Int, j: Int, k: Int): Boolean {
        return (potential[k][j] + data.travel_times[0][graph[j]]) <
                (potential[k][i] + data.travel_times[0][graph[i]] +
                sumDistance[j+1] - sumDistance[i+1] + 0.0001F)
    }

}



class Individual(private val data: Data,
                 private val initGene: List<Int>? = null) {

    private val numNurses = data.nbr_nurses
    private val numPatients = data.patients.size


    var gene: List<Int> = initGene ?: (IntArray(numPatients){ it } + IntArray(numNurses-1){ -(it+1) }).asList().shuffled()
    fun saveToFile() {
        val solutionGene = mutableListOf<Int>()
        gene.forEach { solutionGene.add(it+1) }
        val solution = splitListOnNegative(solutionGene, 1)
        println(solution)
    }
    fun getFitness(penalty: Float = 0F, route: List<Int> = gene): Float {
        return getTravelDistance(route) + penalty * getTimeWindowViolation(route) + getCapacityViolation(route)
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
        var travelTime: Float
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

                // is total time used less then end time
                timeViolation += max(0F, usedCareTime - data.patients[newLocation].end_time.toFloat())
                prevLocation = newLocation
            }
        }

        return timeViolation + max(0F, timeViolation - data.depot.return_time)
    }
    private fun getCapacityViolation(route: List<Int> = gene): Float {
        var sumDemamd = 0F
        var capacitySum = 0F
        var prevLocation = 0
        for (patient in route) {
            // If route is finished or not
            if (patient < 0) {
                sumDemamd = 0F
            }
            else {
                sumDemamd += data.patients[patient].demand
                if (sumDemamd > data.capacity_nurse)
                    capacitySum += (sumDemamd - data.capacity_nurse)
            }

        }
        return capacitySum
    }

    fun mutate(mutateProbability: Float = 0.3F) {
        // chance of mutating
        if (Random.nextFloat() < mutateProbability) {
            when (Random.nextInt(0,3)) {
                0 -> swapRoutePatients()
                1 -> insertMutation()
                2 -> swapTwoBetweenRoutes()
            }
        }
    }
    fun movePatientToRoute() {
        val routes = splitListOnNegative(gene)
        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var bestFit = Pair(listOf(-1,-1, listOf<Int>()), 0F)

        println("\nRoute1 $route")
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

            if (neighborRoute.isEmpty()) error("Neighbor route is EMPTY")
            val neighborIndex = neighborRoute.indexOf(shortest.first)
            bestFit = Pair(listOf(-1,-1, listOf<Int>()), getFitness(10F, neighborRoute))

            println("Neighbor route $neighborRoute")
            // checks the 2± neighborhood of the best neighboring route node
            for (i in maxOf(0, neighborIndex-2) .. minOf(neighborIndex+2, neighborRoute.size)) {

                val tempRoute = neighborRoute.toMutableList()
                tempRoute.add(i, node)
                val fitness = getFitness(10F, tempRoute)
                if (fitness < bestFit.second) {
                    println("New best fitness")
                    bestFit = Pair(listOf(node, i, neighborRoute), fitness)
                }


            }
            if (bestFit.first[0] != -1)
                break
        }
        if (bestFit.first[0] != -1) {
            val (node, index, inRoute) = bestFit.first
            if (node !is Int) error("node is not an int")
            if (index !is Int) error("index is not an int")
            route.toMutableList().remove(node)
            routes.add(route)

            val newRoute = inRoute as MutableList<Int>
            routes.remove(newRoute)
            print("Switched parents. Old fitness:${getFitness(10F, newRoute)}")
            newRoute.add(index, node)
            routes.add(newRoute)
            println(", New fitness:${getFitness(10F, newRoute)}")

        }
        else {
            routes.add(route)
            println("No nodes switched between routes")
        }
        concatenateGene(routes)

    }

    private fun swapTwoBetweenRoutes() {

        val routes = splitListOnNegative(gene)
        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var bestFit = Pair(listOf(-1,-1, listOf<Int>()), getFitness(10F, route))


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

            if (neighborRoute.isEmpty()) error("Neighbor route is EMPTY")

            //println("Route2 $neighborRoute")

            var bestNeightborFitness = getFitness(10F, neighborRoute)

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

                val fitness1 = getFitness(10F, temproute1)
                val fitness2 = getFitness(10F, temproute2)
                // better for route1 and better for route2
                if (fitness1 < bestFit.second && fitness2 < bestNeightborFitness) {
                    bestNeightborFitness = fitness2
                    bestFit = Pair(listOf(node,neighborNode, neighborRoute), fitness1)
                }
            }

            if (bestFit.first[0] != -1) {
                println("Init Route1 fitness ${getFitness(10F, route)}, Route2 fitness: ${getFitness(10F,neighborRoute)}")
                println("Best fit, Node ${bestFit.first[0]} and ${bestFit.first[1]}")
            }

        }

        if (bestFit.first[0] != -1) {
            println("\nBest overall: ${bestFit.first[0]} with ${bestFit.first[1]}, fitness ${bestFit.second}")

            val newRoute1 = route.toMutableList()
            val newRoute2 = bestFit.first[2] as MutableList<Int>
            routes.remove(newRoute2)
            val node1 = newRoute1.indexOf(bestFit.first[0])
            val node2 = newRoute2.indexOf(bestFit.first[1])
            newRoute1[node2] = bestFit.first[0] as Int
            newRoute2[node1] = bestFit.first[1] as Int
            routes.add(newRoute1)
            routes.add(newRoute2)

        }
        else {
            routes.add(route)
            print("No good switch")
        }

        // route to gene and save it
        concatenateGene(routes)

    }
    private fun swapRoutePatients() {
        val routes = splitListOnNegative(gene)

        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        var reducedRoute = route.toMutableList()
        var bestFit = Pair(listOf(-1,-1), Float.MAX_VALUE)

        for (nodePos1 in route.indices) {
            for (nodePos2 in (nodePos1+1) until route.size) {
                Collections.swap(reducedRoute, nodePos1, nodePos2)
                val fitness = getFitness(10F, reducedRoute) // TODO("Make time penalty available for Individual class")
                if (fitness < bestFit.second)
                    bestFit = Pair(listOf(nodePos1, nodePos2), fitness)
                reducedRoute = route.toMutableList() // back to original
            }
            //if (bestFit.first[0] != -1)
            //    println("Best swap, ${bestFit.first[0]} with ${bestFit.first[1]}, fitness ${bestFit.second}")
        }

        if (bestFit.first[0] != -1)
            Collections.swap(reducedRoute, bestFit.first[0], bestFit.first[1])

        if (bestFit.second < getFitness(10F, route))
            routes.add(reducedRoute)
        else
            routes.add(route)

        println("Original route $route, Fitness: ${getFitness(10F, route)}")
        println("New route      $reducedRoute, Fitness: ${getFitness(10F, reducedRoute)}")

        // route to gene and save it
        concatenateGene(routes)
    }
    private fun insertMutation() {

        val routes = splitListOnNegative(gene)

        var route: List<Int>
        do {
            route = routes.random()
        } while (route.isEmpty())
        routes.remove(route)

        val reducedRoute = route.toMutableList()
        for (node in route) {
            val nodePosition = route.indexOf(node)

            reducedRoute.remove(node)
            var bestFit: Pair<Int, Float> = Pair(-1, Float.MAX_VALUE)

            for (position in 0 .. reducedRoute.size) {
                val newRoute = reducedRoute.toMutableList() // make copy of list
                newRoute.add(position, node)
                val fitness = getFitness(10F, newRoute)
                //println(fitness)
                if (fitness < bestFit.second)
                    bestFit = Pair(position, fitness)

            }

            // if no feasible position: insert node back into original spot
            if (bestFit.first == -1)
                reducedRoute.add(nodePosition, node)
            else
                reducedRoute.add(bestFit.first, node)

        }
        println("Original route $route, Fitness: ${getFitness(10F, route)}")
        println("New route      $reducedRoute, Fitness: ${getFitness(10F, reducedRoute)}")
        routes.add(reducedRoute)
        // route to gene and save it
        concatenateGene(routes)
    }
    private fun swapMutation() {
        val point1 = Random.nextInt(0,gene.size)
        var point2 = Random.nextInt(0,gene.size)
        while (point1 == point2) point2 = Random.nextInt(0,gene.size)
        Collections.swap(gene, point1, point2)
    }
    private fun inversionMutation() {
        val point1 = Random.nextInt(0,gene.size)
        var point2 = Random.nextInt(0,gene.size)
        while (point1 == point2) point2 = Random.nextInt(0,gene.size)
        if (point2 < point1)
            gene.toMutableList().subList(point2,point1).reverse()
        else
            gene.toMutableList().subList(point1,point2).reverse()

        gene.toList()
    }

    fun removeCluster() {
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

            println("$previousNode to $nextNode, distance $arcDistance")

            if (arcDistance > longestArc) {
                longestArc = arcDistance
                longestArcIndex = nextNodeIndex
            }
            previousNode = nextNode
        }
        println("$selectedRoute, Longest arc Index: $longestArcIndex")

        val shortendRoute = selectedRoute.subList(0, longestArcIndex)
        val unassignedNodes = selectedRoute.subList(longestArcIndex, selectedRoute.lastIndex+1)
        // add shortened route back in
        routes.add(shortendRoute)

        while (unassignedNodes.isNotEmpty()) {


        }

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
        //println("BEfore loop")
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
    /*fun getEdges(): MutableMap<Int, MutableList<Int>> {
        val edges = mutableMapOf<Int, MutableList<Int>>()

        val graph = mutableListOf<Int>()

        // make graph and store depot indices
        for (i in gene.indices) {
            if (gene[i] >= 0) {
                graph.add(gene[i])
            }
            else {
                splitIndexes.add(i) // indexes of depot
            }
        }

        for (pos in graph.indices) {
            edges[graph[pos]] = when (pos) {
                0 -> mutableListOf<Int>(graph.last(), graph[1])
                (graph.size - 1) -> mutableListOf<Int>(graph[pos - 1], graph[1])
                else -> {
                    mutableListOf<Int>(graph[pos - 1], graph[pos + 1])
                }
            }
        }
        return edges
    }*/
    fun addDelimiters() {
        // make graph and store depot indices
        //val split = LinearBoundedSplit(data, gene as MutableList<Int>)
        //gene = split.findDelimiters()
        print("")
    }

    fun insertionHeuristicSimple(used: List<Int>) {

        // patientStartTime = [ (index, patient), ... ], index = patientID aka node number
        val patientList = data.patients.mapIndexed { index: Int, patient -> index to patient  }.toMutableList()

        var usedNurses = 1
        val initRoutes = mutableListOf<List<Int>>()
        var firstRun = true

        while (patientList.isNotEmpty() && (usedNurses < (data.nbr_nurses-2) )) {

            var nodePair = patientList.sortedBy { it.second.start_time }.toMutableList().first()

            if (firstRun) {
                var i = 0
                while (nodePair.first in used) {
                    nodePair = patientList.sortedBy { it.second.start_time }.toMutableList()[i++]
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
        concatenateGene(initRoutes)
    }
    fun insertionHeuristic(used: List<Int>) {
        // TODO("Fix so that each route does not end on a bad node")
        // patientStartTime = [ (index, patient), ... ], index = patientID aka node number
        val patientList = data.patients.mapIndexed { index: Int, patient -> index to patient  }.toMutableList()


        var usedNurses = 1
        val initRoutes = mutableListOf<List<Int>>()
        var firstRun = true


        while (patientList.isNotEmpty() && (usedNurses < (data.nbr_nurses-2) )) {

            var nodePair = patientList.sortedBy { it.second.start_time + it.second.care_time }.toMutableList().first()

            if (firstRun) {
                var i = 0
                while (nodePair.first in used) {
                    nodePair = patientList.sortedBy { it.second.start_time + it.second.care_time }.toMutableList()[i++]
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
                            val fitness =  getTravelDistance(temproute) + node.second.start_time
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
            initRoutes.add(route)
        }
        // add remaining un-assigned patients
        val route = mutableListOf<Int>()
        patientList.forEach { route.add(it.first) }
        initRoutes.add(route)

        concatenateGene(initRoutes) // save gene
    }
    private fun splitListOnNegative(arr: List<Int>, onValue: Int = 0): MutableList<List<Int>> {
        // splits on values less than 0
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
    private fun concatenateGene(routes: List<List<Int>>) {
        val newGene = mutableListOf<Int>()

        var nurse = -1
        for (route in routes) {
            if (route.isNotEmpty()) {
                route.forEach { newGene.add(it)}
                newGene.add(nurse--)
            }
        }
        // adds remaining nurses to gene
        val nursesLeft = data.nbr_nurses + nurse
        repeat(nursesLeft) {
            newGene.add(nurse--)
        }
        gene = newGene
    }
}



data class GenerationData(var ptype: MutableList<Float> = mutableListOf<Float>(),
                          var fittest: MutableList<Individual> = mutableListOf<Individual>(),
                          var meanFitness: MutableList<Double> = mutableListOf())


class GeneticAlgorithm(private val data: Data,
                       private val sizeOfPopulation: Int) {

    var population = mutableListOf<Individual>()
    private val tournamentSize = 2
    private val mutateProbability = 0.8F
    private val timeWindowPenalty = 5F
    var dataStorage = GenerationData()
    lateinit var fittestIndividual: Individual

    private fun storeCurrentState() {
        dataStorage.ptype.add(uniquePhenoTypes())
        //dataStorage.meanFitness.add(getAverageFitness())
    }

    private fun uniquePhenoTypes(): Float {
        val fitness = mutableListOf<Float>()
        population.forEach { fitness.add(it.getFitness()) }

        return fitness.distinct().size / fitness.size.toFloat()
    }

    private fun initPopulation() {

        // repeat(sizeOfPopulation) { population.add(Individual(this.data)) }
        //kmeansInit()
        constructionHeuristic()
    }
    fun kmeansInit() {
        val xyPosition = getPatientCoords()
        //val clusters = xmeans(patientPositions, data.nbr_nurses)
        val xPos = normalizeData(xyPosition[0])
        val yPos = normalizeData(xyPosition[1])
        //println("Number of clusters ${clusters.k}")
        //val kmeansK = kmeans(patientPositions, data.nbr_nurses).k
        val xyCoords = Array<DoubleArray>(data.patients.size) { doubleArrayOf(0.0) }
        for (i in xPos.indices)
            xyCoords[i] = doubleArrayOf(xPos[i], yPos[i])


        for (i in 2.. data.nbr_nurses) {
            val labels = kmeans(data = xyCoords, k = i)

            //println("Kmeans: Number of clusters ${labels.k}")
            //labels.y.forEach { println(it) }

            val patients = IntArray(data.patients.size) { it }
            val patientMap = patients.groupBy { labels.y[it] }
            //println(patientMap)
            val newGen = mutableListOf<Int>()
            var nurse = -1
            for (route in patientMap.values) {
                route.forEach { newGen.add(it) }
                newGen.add(nurse--)
            }
            population.add(Individual(data, newGen))
        }
        /*
       val newgene = mutableListOf<Int>()
       var routenr = -1
       for (route in patientMap.values)  {
           println(route)
           val dist = mutableListOf<Float>()
           for (i in 1 until  route.size-1) {
               dist.add(data.travel_times[])
           }
           newgene.add(routenr--)
       }*/
    }

    fun constructionHeuristic() {
        var used = mutableListOf<Int>()

        if (true) {
            repeat(sizeOfPopulation) {
                val indiv = Individual(this.data)
                indiv.insertionHeuristic(used)
                used.add(indiv.gene[0])
                population.add(indiv)
            }
        }
        else {
            used = mutableListOf<Int>()
            repeat(sizeOfPopulation) {
                val indiv = Individual(this.data)
                indiv.insertionHeuristicSimple(used)
                used.add(indiv.gene[0])
                population.add(indiv)
            }
        }

        /*
        repeat(sizeOfPopulation-15) {
            population.add(Individual(this.data))
        }
        */

    }

    private fun normalizeData(data: Array<Double>): List<Double> {
        val maxValue = data.maxOrNull() ?: error("Must be error in coords")
        val minValue = data.minOrNull() ?: error("Must be error in coords")

        val normalizedData = mutableListOf<Double>()
        data.forEach { normalizedData.add((it - minValue)/(maxValue-minValue)) }
        return normalizedData.toList()


    }
    private fun getPatientCoords(): List<Array<Double>> {

        val xCoord = Array<Double>(data.patients.size){0.0}
        val yCoord = Array<Double>(data.patients.size){0.0}
        for (i in data.patients.indices) {
            xCoord[i] = data.patients[i].x_coord.toDouble()
            yCoord[i] = data.patients[i].y_coord.toDouble()
        }
            //xyCoords[i] = doubleArrayOf(data.patients[i].x_coord.toDouble(), data.patients[i].y_coord.toDouble())

        return listOf(xCoord, yCoord)
    }

    private fun getAverageFitness(): Double {
        val fitness = mutableListOf<Float>()
        population.forEach { fitness.add(it.getFitness(timeWindowPenalty)) }
        val fittestIndex = fitness.indices.minByOrNull{ fitness[it] } ?: error("This should not run")
        fittestIndividual = population[fittestIndex]
        //dataStorage.fittest.add(population[fittestIndex])
        return fitness.average()
    }

    fun  fit(numGenerations: Int) {
        this.initPopulation()
        println("Population of ${population.size} initiated")

        val printNr = if (numGenerations < 30) numGenerations else (numGenerations * 0.01 + 1).roundToInt()
        for (i in 1..numGenerations) {
            doGeneration()

            if (i % printNr == 0) {
                //storeCurrentState()
                println("Gen: $i, avgFitness: ${getAverageFitness()}, Fittest: ${fittestIndividual.getFitness()}")
            }
        }

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

    }

    private fun parentSelection(): Array<Individual> {
        val mother: Individual = getParent()
        var father: Individual = getParent()
        while (mother == father) {  father = getParent() }
        return arrayOf(mother, father)
    }

    private fun getParent(): Individual {
        // TODO("Add more advanced parent selection methods")
        return tournamentSelection()
    }

    private fun tournamentSelection(): Individual {

        var bestParent = population.random()
        var competingParent = population.random()
        repeat(tournamentSize) {
            while (bestParent == competingParent) { competingParent = population.random() }
            if (bestParent.getFitness(timeWindowPenalty) > competingParent.getFitness(timeWindowPenalty)) { // minimizing fitness
                bestParent = competingParent
            }
        }
        return bestParent
    }

    private fun getOffspring(mother: Individual, father: Individual): List<Individual> {
        val (firstChild, secondChild) = doCrossover(mother, father)
        firstChild.mutate(mutateProbability)
        secondChild.mutate(mutateProbability)
        return listOf(firstChild, secondChild)

    }

    private fun doCrossover(mother: Individual, father: Individual): List<Individual> {
        // TODO("Add more advanced crossover methods")

        //if (Random.nextFloat() < 1) {
        //    mother.removeCluster()
        //    father.removeCluster()
        //    return listOf(mother, father)
        //}

        val firstChild = alternatingEdgesCrossover(mother, father)
        val secondChild = alternatingEdgesCrossover(father, mother)
        return listOf(firstChild, secondChild)

    }

    private fun toGraph(gene: List<Int>): List<Int> {
        val graph = mutableListOf<Int>()
        gene.forEach{if (it >= 0) {graph.add(it)} }
        return graph
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

            if (nextArc in childGene) {
                nextArc = unusedEdges.random()
            }
            childGene.add(nextArc)
            unusedEdges.remove(nextArc)
            if (unusedEdges.isEmpty()) {break}
        }
        //val newGene = simpleSplitToRoutes(childGene)
        //val child = Individual(data, childGene)
        //child.addDelimiters() // adds de-limiters (nurses)

        return Individual(data, childGene)
    }

    private fun simpleSplitToRoutes(graph: List<Int>): List<Int> {
        var num_vehicle = 0
        var remained_cap = data.capacity_nurse //- 100
        val newGene: MutableList<Int> = mutableListOf<Int>()
        for (node in graph) {
            if ((remained_cap - data.patients[node].demand) >= 0) {
                newGene.add(node)
                remained_cap -= data.patients[node].demand
            }
            else {
                newGene.add( -(num_vehicle+1) )
                newGene.add(node)
                num_vehicle += 1
                remained_cap = data.capacity_nurse - data.patients[node].demand
            }
        }

        if ((graph.size + (data.nbr_nurses-1) - newGene.size) != 0) {
            newGene.addAll(MutableList((graph.size + data.nbr_nurses-1) - newGene.size){ -(it+num_vehicle+1) })
        }


        return newGene

    }
    private fun oneOrderCrossover(mother: Individual, father: Individual): Individual {
        // TODO("Do crossover")
        val rndNumbers = listOf<Int>(Random.nextInt(0, mother.gene.size), Random.nextInt(0, mother.gene.size)).sorted()
        val firstChild = mother.gene.toMutableList() // simple copy


        //firstChild = mother.gene.subList(rndNumbers[0], rndNumbers[1])
        return Individual(data, firstChild)
    }

    private fun survivorSelection(offsprings: List<Individual>): List<Individual> {
        // TODO("Add more complex survivor selection methods")
        return nFittest(offsprings)
    }

    private fun nFittest(offsprings: List<Individual>): List<Individual> {
        // TODO("Add elitism")
        val competingPopulation = offsprings + population

        val fitnessList = competingPopulation.associateBy({it.getFitness()}, {it})
        val fittestIndividuals = fitnessList.toSortedMap().values.toList().subList(0, population.size)


        return fittestIndividuals
    }
}





fun main(args: Array<String>) {

    val filename = "train_0.json"
    val data = createDataclass(filename)

    val model = GeneticAlgorithm(data, 4)

    val testing = true
    if (!testing) {
        val timeInMillis = measureTimeMillis {
            model.fit(10)
        }

        println("(The operation took $timeInMillis ms)")
        println("\nThe fittest individual: ${model.fittestIndividual.getFitness(0F)}")
        println("\nThe fittest individual with constraint: ${model.fittestIndividual.getFitness(1F)}")
        model.fittestIndividual.saveToFile()
    }

    model.constructionHeuristic()
    model.population[0].gene = listOf(19, 23, 24, 26, 28, 29, 27, 25, 22, 21, 20, 74, -1, 66, 64, 62, 61, 73, 71, 60, 63, 67, 65, 68, -2, 4, 2, 6, 7, 9, 10, 8, 5, 3, 1, 0, 46, -3, 42, 41, 40, 39, 58, 43, 45, 44, 47, 50, 49, 51, 48, -4, 89, 86, 85, 82, 81, 83, 84, 87, 88, 90, -5, 12, 16, 17, 18, 14, 15, 13, 11, 98, -6, 97, 95, 94, 93, 91, 92, 96, 99, -7, 31, 32, 30, 34, 36, 37, 38, 35, 33, -8, 56, 54, 53, 52, 55, 57, 59,  -9, 80, 77, 75, 70, 69, 72, 76, 78, 79, -10, -11, -12, -13, -14, -15, -16, -17, -18, -19, -20, -21, -22, -23, -24)
    repeat(50) {
        model.population[0].movePatientToRoute()
    }

    //model.population.forEach { it.movePatientToRoute() }
    //val testGene = listOf<Int>(4, 2, 6, 9, 7, 10, 8, 5, 3, 1, -1, 12 ,13 ,14 ,15 ,-2 ,-3)
    //Individual(data, testGene).swapTwoBetweenRoutes()
    model.population[0].saveToFile()





}
