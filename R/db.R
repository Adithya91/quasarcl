#' @import rmongodb

#' @export
getDbConnection<-function(host = "localhost" , username = "" , password = "" , db = DATABASE) 
{
	dbConnection <- mongo.create(host, username, password, db)
	
	if (!mongo.is.connected(dbConnection)) 
	{
		stop("MongoDB connection not established")
	}
	return (dbConnection)
}



#' @export
closeDbConnection <- function(connection) 
{
	mongo.destroy(connection)
}



#' @export
saveQuasarSet<-function(connection, quasarSetName, quasarSetDate, db = DATABASE) 
{
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")
	query<-buildQuasarSetData(quasarSetName, quasarSetDate)
	mongo.insert(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



#' @export
getQuasarSetByName<-function(connection, name, db = DATABASE) 
{
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")
	query<-mongo.bson.from.list(list(name = name))
	quasarSet<-mongo.find.one(connection, ns, query)

	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
	if (is.null(quasarSet))
	{
		stop("QuasarSet with such name does not exist")
	}
	return (mongo.bson.to.list(quasarSet))
}



#' @export
getAllQuasarSets<-function(connection, db = DATABASE) 
{
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")
	quasarSets<-mongo.cursor.to.list(mongo.find(connection, ns))
	return (quasarSets)
}



#' @export
deleteQuasarSet<-function(connection, quasarSetName, db = DATABASE)
{
	#usuwanie kwazarów z danego zbioru
	deleteQuasars(connection, quasarSetName, db)	
	
	#usuwanie zbioru
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName, db)
	query<-mongo.bson.from.list(list('_id' = quasarSet$'_id'))
	mongo.remove(connection, ns, query)	
}



#' @export
saveQuasars<-function(connection, quasars, quasarSetName, db = DATABASE ) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName, db)
	for(quasar in quasars) 
	{
		query<-buildQuasarData(quasar, quasarSet$'_id')
		mongo.insert(connection, ns, query)
		
		err <- mongo.get.last.err(connection, db)
		if (!is.null(err)) 
		{
			stop(mongo.get.server.err.string(connection))
		}
	}
	updateQuasarSet(connection, quasarSet$'_id', length(quasars), db)
}



#' @export
getQuasars<-function(connection, quasarSetName, db = DATABASE, n = 0L) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName, db)
	query<-mongo.bson.from.list(list(quasarSetOID = quasarSet$'_id'))
	cursor<-mongo.find(connection, ns, query, limit = n)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err))
	{
		stop(mongo.get.server.err.string(connection))
	}
	return (mongo.cursor.to.list(cursor))
}



#' @export
saveParamResultSet<-function(connection, name, date, options, db = DATABASE) 
{
	ns<- paste(DATABASE, PARAM_RESULTSET_COLLECTION, sep=".")
	query<-buildParamResultSetData(name, date, options)
	mongo.insert(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}


#' @export
getParamResultSetByName <- function(connection, resultSetName, db = DATABASE) 
{	
	ns<- paste(db, PARAM_RESULTSET_COLLECTION, sep=".")
	query<-mongo.bson.from.list(list(name = resultSetName))
	resultSet<-mongo.find.one(connection, ns, query)

	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
	if (is.null(resultSet))
	{
		stop("QuasarSet with such name does not exist")
	}
	return (mongo.bson.to.list(resultSet))
}



#' @export
deleteParamResultSet <- function(connection, resultSetName, db = DATABASE)
{
	#usuwanie wyników z danego zbioru
	deleteParamResults(connection, resultSetName, db)
	
	#usuwanie zbioru
	ns<- paste(db, PARAM_RESULTSET_COLLECTION, sep=".")
	resultSet <- getParamResultSetByName(connection, resultSetName, db)
	query<-mongo.bson.from.list(list('_id' = resultSet$'_id'))
	mongo.remove(connection, ns, query)
}



#' @export
saveParamResults<-function(connection, resultSetName, results, db = DATABASE) 
{
	ns<- paste(db, PARAM_RESULT_COLLECTION, sep=".")
	resultSet <- getParamResultSetByName(connection, resultSetName, db)
	query <- lapply(results, buildParamResultData, resultSet$'_id')
	mongo.insert.batch(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
	updateParamResultSet(connection, resultSet$'_id', length(results), db)
}


#' @export
getParamResults <- function(connection, resultSetName, db = DATABASE)
{
	ns<- paste(db, PARAM_RESULT_COLLECTION, sep=".")
	resultSet <- getParamResultSetByName(connection, resultSetName, db)
	query<-mongo.bson.from.list(list(paramResultSetOID = resultSet$'_id'))
	cursor<-mongo.find(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err))
	{
		stop(mongo.get.server.err.string(connection))
	}
	return (mongo.cursor.to.list(cursor))
}



updateQuasarSet<-function(connection, oid, size, db = DATABASE) 
{
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")

	buffer <- mongo.bson.buffer.create()
	mongo.bson.buffer.start.object(buffer, "$inc")
	mongo.bson.buffer.append.long(buffer, "size", size)
	mongo.bson.buffer.finish.object(buffer)
	query <- mongo.bson.from.buffer(buffer)
	
	mongo.update(connection, ns, list('_id' = oid), query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



deleteQuasars<-function(connection, quasarSetName, db = DATABASE) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName, db)
	query<-mongo.bson.from.list(list(quasarSetOID = quasarSet$'_id'))
	mongo.remove(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



updateParamResultSet<-function(connection, oid, size, db = DATABASE) 
{
	ns<- paste(db, PARAM_RESULTSET_COLLECTION, sep=".")

	buffer <- mongo.bson.buffer.create()
	mongo.bson.buffer.start.object(buffer, "$inc")
	mongo.bson.buffer.append.long(buffer, "size", size)
	mongo.bson.buffer.finish.object(buffer)
	query <- mongo.bson.from.buffer(buffer)
	
	mongo.update(connection, ns, list('_id' = oid), query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



deleteParamResults <- function(connection, resultSetName, db = DATABASE)
{
	ns<- paste(db, PARAM_RESULT_COLLECTION, sep=".")
	resultSet <- getParamResultSetByName(connection, resultSetName, db)
	query<-mongo.bson.from.list(list(paramResultSetOID = resultSet$'_id'))
	mongo.remove(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



buildQuasarSetData<-function(quasarSetName, quasarSetDate)
{
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "name", quasarSetName)
	mongo.bson.buffer.append(buffer, "date", quasarSetDate)
	return (mongo.bson.from.buffer(buffer))	
}



buildQuasarData <- function(quasar, quasarSetID) 
{
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "quasarSetOID", quasarSetID)
	mongo.bson.buffer.append(buffer, "params", quasar$params)
	mongo.bson.buffer.append(buffer, "size", length(quasar$values))
	mongo.bson.buffer.append(buffer, "values", quasar$values)
	mongo.bson.buffer.append(buffer, "error", quasar$error)
	return (mongo.bson.from.buffer(buffer))
}



buildParamResultSetData<-function(name, date, options)
{
	names(options$continuumWindows) <- NULL
	names(options$feWindows) <- NULL
	
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "name", name)
	mongo.bson.buffer.append(buffer, "setName", paste(name, as.numeric(date), sep=""))
	mongo.bson.buffer.append(buffer, "date", date)
	mongo.bson.buffer.append(buffer, "options", options)
	return (mongo.bson.from.buffer(buffer))
}



buildParamResultData <- function(result, resultSetID) 
{
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "paramResultSetOID", resultSetID)
	mongo.bson.buffer.append(buffer, "mjd", result$mjd)
	mongo.bson.buffer.append(buffer, "plate", result$plate)
	mongo.bson.buffer.append(buffer, "fiber", result$fiber)
	mongo.bson.buffer.append(buffer, "continuumChisq", result$continuumChisq)
	mongo.bson.buffer.append(buffer, "continuumReglin", result$continuumReglin)
	mongo.bson.buffer.append(buffer, "reglin", result$reglin)
	mongo.bson.buffer.append(buffer, "feScaleRate", result$feScaleRate)
	mongo.bson.buffer.append(buffer, "feWindowsSize", result$feWindowsSize)
	mongo.bson.buffer.append(buffer, "feWindowsReducedChisq", result$feWindowsReducedChisq)
	mongo.bson.buffer.append(buffer, "feFullReducedChisq", result$feFullReducedChisq)
	mongo.bson.buffer.append(buffer, "feFullEW", result$feFullEW)
	mongo.bson.buffer.append(buffer, "feRangeReducedChisq", result$feRangeReducedChisq)
	mongo.bson.buffer.append(buffer, "feRangeEW", result$feRangeEW)
	mongo.bson.buffer.append(buffer, "elementsFits", result$elementsFits)
	
	return (mongo.bson.from.buffer(buffer))
}
