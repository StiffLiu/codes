DELIMITER $$
DROP PROCEDURE find_scc $$
CREATE PROCEDURE find_scc() 
COMMENT 'Find Strongly Connected component from table scc'
BEGIN
   DECLARE var_vertex_stack_count INT DEFAULT 0;
   DECLARE var_component INT DEFAULT 0;
   DECLARE var_call_stack_index INT DEFAULT 0;
   DECLARE var_vertex_index INT DEFAULT 0;
   DECLARE var_vertex_index_global INT DEFAULT 0;
   -- vertex from which we start the construction process
   DECLARE var_vertex INT DEFAULT NULL;
   DECLARE var_adj INT DEFAULT NULL;
   DECLARe var_ret_value INT DEFAULT NULL;
   
   DECLARE var_adj_link_num INT DEFAULT NULL;
   DECLARE var_adj_in_vertex_stack INT DEFAULT NULL;
   
   DECLARE var_helper INT DEFAULT 0;
   DECLARE var_helper1 INT DEFAULT 1;
   
   -- CREATE TEMPORARY TABLE temp_graph (src int not null, dest int not null, primary key (src, dest));
   -- CREATE TEMPORARY TABLE temp_vertex_stack (idx int primary key, vertex int not null));
   -- CREATE TEMPORARY TABLE temp_link (vertex int primary key, num int, in_vertex_stack int));
   -- CREATE TEMPORARY TABLE temp_call_stack (vertex int not null, adj int not null, idx int primary key, vertex_index int  not null));
   -- CREATE TEMPORARY TABLE temp_scc (vertex int primary key, component int not null));
   -- CREATE TABLE src_graph (src int not null, dest int not null, primary key(src, dest)));

   DELETE FROM temp_graph;
   DELETE FROM temp_vertex_stack;
   DELETE FROM temp_scc;
   DELETE FROM temp_call_stack;
   INSERT temp_graph(src, dest) SELECT distinct src, dest FROM src_graph where src != dest;
   
   -- try to remove more edges to avoid unnecessary checks.
   SELECT COUNT(*) INTO var_helper1 FROM temp_graph;
   WHILE var_helper != var_helper1 DO
		-- SELECT * FROM temp_graph;
		SET var_helper = var_helper1;
		DELETE FROM temp_link;
		INSERT INTO temp_link(vertex, num, in_vertex_stack) SELECT DISTINCT dest, NULL, NULL FROM temp_graph WHERE dest IN (SELECT DISTINCT src FROM temp_graph);
		-- ignore vertices whose in-degree or out-degree is zero, which cannot appear in any loop.
		DELETE FROM temp_graph WHERE (src NOT IN (SELECT DISTINCT vertex FROM temp_link) ) OR (dest NOT IN (SELECT DISTINCT vertex FROM temp_link) ) ;
		SELECT COUNT(*) INTO var_helper1 FROM temp_graph;
   END WHILE;
   
   -- make a back up of temp_graph here for efficiently recover the edges in cycles.
   
   DELETE from temp_graph where src = dest;
   
   SELECT src INTO var_vertex FROM temp_graph LIMIT 0,1;
   WHILE var_vertex IS NOT NULL DO
		SET var_adj = NULL;
		SET var_call_stack_index = 0;
		SET var_vertex_stack_count = 0;
		-- SELECT "Begin vertex", var_vertex;
		-- SELECT * FROM temp_graph;
	begin_of_find:
		WHILE var_call_stack_index >=0 DO
			IF var_adj IS NULL THEN
				-- first call;
				UPDATE temp_link SET num = var_vertex_index_global, in_vertex_stack = 1 WHERE vertex = var_vertex;
				INSERT temp_vertex_stack(idx, vertex) VALUES(var_vertex_stack_count, var_vertex);
				SET var_vertex_index = var_vertex_index_global;
				SET var_vertex_index_global = var_vertex_index_global + 1;
				SET var_vertex_stack_count = var_vertex_stack_count + 1;
			ELSE
				-- return from previous call
				UPDATE temp_link SET num = var_ret_value WHERE vertex = var_vertex AND var_ret_value < num;
				SET var_adj = NULL;
			END IF;
			
			-- adj is NULL at this point, pick the next available adjacency vertex.
			SELECT dest INTO var_adj FROM temp_graph WHERE src = var_vertex LIMIT 0,1;
			-- SELECT var_vertex, var_vertex_index, var_adj;
			WHILE var_adj IS NOT NULL DO
				-- if there exists one adjacency vertex, visit it.
				-- delete the edge to avoid being revisited.				
				DELETE FROM temp_graph WHERE src = var_vertex AND dest = var_adj;
				SET var_adj_link_num = NULL;
				SET var_adj_in_vertex_stack = NULL;
				SELECT num, in_vertex_stack INTO var_adj_link_num, var_adj_in_vertex_stack FROM temp_link WHERE vertex = var_adj;
				-- SELECT var_adj_link_num, var_adj_in_vertex_stack;
				IF var_adj_link_num IS NULL THEN
					-- if this adjacency vertex hasn't already visited
					-- save current to call stack
					INSERT INTO temp_call_stack(idx, vertex, adj, vertex_index) VALUES(var_call_stack_index, var_vertex, var_adj, var_vertex_index);
					SET var_call_stack_index = var_call_stack_index + 1;
					SET var_vertex = var_adj;
					SET var_adj = NULL;
					ITERATE begin_of_find;
				ELSEIF var_adj_in_vertex_stack IS NOT NULL THEN
					UPDATE temp_link SET num = var_adj_link_num WHERE vertex = var_vertex AND var_adj_link_num < num;
				END IF;
				SET var_adj = NULL;
				SELECT dest INTO var_adj FROM temp_graph WHERE src = var_vertex LIMIT 0,1;
				-- SELECT var_vertex, var_vertex_index, var_adj;
			END WHILE;			
			
							
			-- all the adjacency vertices has already been visited;
				
			-- check if the lowest link of the vertex is the same as that of it's 
			-- vertex_index, if yes a connected component has been found.
			SELECT COUNT(*) INTO var_helper FROM temp_link WHERE vertex = var_vertex AND num = var_vertex_index;
			IF var_helper > 0 THEN
				-- SELECT "Find one component", var_vertex_index;
				-- a strongly connected component has been found.
				SELECT idx INTO var_vertex_stack_count FROM temp_vertex_stack WHERE vertex = var_vertex;
				INSERT INTO temp_scc(vertex, component) SELECT vertex, var_component FROM temp_vertex_stack WHERE idx >= var_vertex_stack_count;
				UPDATE temp_link SET in_vertex_stack = NULL WHERE vertex IN (SELECT DISTINCT vertex FROM temp_vertex_stack WHERE idx >= var_vertex_stack_count);
				DELETE FROM temp_vertex_stack WHERE idx >= var_vertex_stack_count;
				SET var_component = var_component + 1;
			END IF;
			
			-- this call is finished, retun value in ret_value
			-- and recover call stack.
			SELECT num INTO var_ret_value FROM temp_link WHERE vertex = var_vertex;
			SET var_call_stack_index = var_call_stack_index - 1;
			SELECT vertex, adj, idx, vertex_index INTO var_vertex, var_adj, var_call_stack_index, var_vertex_index FROM temp_call_stack WHERE idx = var_call_stack_index;
			DELETE FROM temp_call_stack WHERE idx = var_call_stack_index;
		END WHILE;
		
		SET var_vertex = NULL;
		SELECT src INTO var_vertex FROM temp_graph LIMIT 0,1;
   END WHILE;
END $$
CALL find_scc $$

DELIMITER ;

