#' @import rmongodb


#' @export
countQuasarsByPeaks <- function(connection, resultSetName, db = DATABASE) {
	ns<- paste(DATABASE, PARAM_RESULT_COLLECTION, sep=".")
	resultSet <- getParamResultSetByName(connection, resultSetName, db)
	print(resultSet$'_id')
	pipeline <- list(
		mongo.bson.from.list(list(
			"$match" = list(
				'paramResultSetOID' = resultSet$"_id")
		)),
		mongo.bson.from.list(list(
			"$unwind" = "$elementsFits"
		)),
		mongo.bson.from.list(list(
			"$group" = list(
				"_id" = "$_id",
				"sum" = list( "$sum" = 1 ))
		)),
		mongo.bson.from.list(list(
			"$group" = list(
				"_id" = "$sum",
				"quasars" = list("$sum" = 1))
		)),
			mongo.bson.from.list(list(
			"$project" = list(
				"_id" = 0,
				"peaks" = "$_id",
				"quasars" = 1)
		))
	)
	result <- mongo.aggregation(connection, ns, pipeline);
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
	return(mongo.bson.to.list(result))
}