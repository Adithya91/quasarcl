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
getQuasars<-function(connection, quasarSetName, db = DATABASE, n = 0L) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName)
	query<-mongo.bson.from.list(list(quasar_set_oid = quasarSet$'_id'))
	cursor<-mongo.find(connection, ns, query, limit = n)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err))
	{
		stop(mongo.get.server.err.string(connection))
	}
	return (mongo.cursor.to.list(cursor))
}



#' @export
deleteQuasars<-function(connection, quasarSetName, db = DATABASE) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName)
	query<-mongo.bson.from.list(list(quasar_set_oid = quasarSet$'_id'))
	mongo.remove(connection, ns, query)
	
	err <- mongo.get.last.err(connection, db)
	if (!is.null(err)) 
	{
		stop(mongo.get.server.err.string(connection))
	}
}



#' @export
saveQuasars<-function(connection, quasars, quasarSetName, db = DATABASE ) 
{
	ns<- paste(db, QUASAR_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName)
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
	updateQuasarSet(connection, quasarSet$'_id', length(quasars))
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
removeQuasarSet<-function(connection, quasarSetName, db = DATABASE)
{
	#usuwanie kwazarów z danego zbioru
	deleteQuasars(connection, quasarSetName)	
	
	#usuwanie zbioru
	ns<- paste(db, QUASARSET_COLLECTION, sep=".")
	quasarSet <- getQuasarSetByName(connection, quasarSetName)
	query<-mongo.bson.from.list(list('_id' = quasarSet$'_id'))
	mongo.remove(connection, ns, query)	
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



buildQuasarSetData<-function(quasarSetName, quasarSetDate)
{
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "name", quasarSetName)
	mongo.bson.buffer.append(buffer, "date", quasarSetDate)
	mongo.bson.buffer.append(buffer, "insertDate", as.POSIXct("2017-03-04"))
	return (mongo.bson.from.buffer(buffer))
	
}



buildQuasarData <- function(quasar, quasarSetID) 
{
	buffer<-mongo.bson.buffer.create()
	oid<-mongo.oid.create()
	mongo.bson.buffer.append(buffer, "_id", oid)
	mongo.bson.buffer.append(buffer, "quasar_set_oid", quasarSetID)
	mongo.bson.buffer.append(buffer, "params", quasar$params)
	mongo.bson.buffer.append(buffer, "size", length(quasar$values))
	mongo.bson.buffer.append(buffer, "values", quasar$values)
	mongo.bson.buffer.append(buffer, "error", quasar$error)
	return (mongo.bson.from.buffer(buffer))
}
