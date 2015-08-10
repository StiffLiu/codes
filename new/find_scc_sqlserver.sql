GO
DROP PROCEDURE find_scc;

GO
-- Find all the edges in strongly connected components.
-- using Tarjan's algorithm.

CREATE PROCEDURE find_scc
AS
BEGIN
   DECLARE @var_vertex_stack_count INT = 0;
   DECLARE @var_component INT = 0;
   DECLARE @var_call_stack_index INT = 0;
   DECLARE @var_vertex_index INT = 0;
   DECLARE @var_vertex_index_global INT = 0;
   -- vertex from which we start the construction process
   DECLARE @var_vertex INT = NULL;
   DECLARE @var_adj INT = NULL;
   DECLARe @var_ret_value INT = NULL;
   
   DECLARE @var_adj_link_num INT = NULL;
   DECLARE @var_adj_in_vertex_stack INT = NULL;
   
   DECLARE @var_helper INT = 0;
   DECLARE @var_helper1 INT = 1;
   
   CREATE TABLE #temp_graph (src INT NOT NULL, dest INT NOT NULL, PRIMARY KEY (src, dest));
   CREATE TABLE #temp_graph_backup (src INT NOT NULL, dest INT NOT NULL, PRIMARY KEY (src, dest));
   CREATE TABLE #temp_vertex_stack (idx INT PRIMARY KEY, vertex INT NOT NULL));
   CREATE TABLE #temp_link (vertex INT PRIMARY KEY, num INT, in_vertex_stack INT));
   CREATE TABLE #temp_call_stack (vertex INT NOT NULL, adj INT NOT NULL, idx int PRIMARY KEY, vertex_index INT NOT NULL));
   CREATE TABLE #temp_scc (vertex INT PRIMARY KEY, component INT NOT NULL));
   -- CREATE TABLE src_graph (src INT NOT NULL, dest INT NOT NULL, PRIMARY KEY(src, dest)));

   DELETE FROM #temp_graph;
   DELETE FROM #temp_vertex_stack;
   DELETE FROM #temp_scc;
   DELETE FROM #temp_call_stack;
   -- change SELECT distinct src, dest FROM src_graph
   -- into your SELECT distcint child_id, parent_id FROM YOUR_TABLE WHERE
   INSERT #temp_graph(src, dest) SELECT DISTINCT src, dest FROM src_graph;
   -- backup self-loops.
   INSERT INTO #temp_graph_backup (src, dest) SELECT src, dest FROM #temp_graph WHERE src = dest;
   DELETE from #temp_graph where src = dest;
   
   -- try to remove more edges to avoid unnecessary checks.
   SELECT @var_helper1 = COUNT(*) FROM #temp_graph;
   WHILE (@var_helper != @var_helper1)
   BEGIN
		-- SELECT * FROM #temp_graph;
		SET @var_helper = @var_helper1;
		DELETE FROM #temp_link;
		INSERT INTO #temp_link(vertex, num, in_vertex_stack) SELECT DISTINCT dest, NULL, NULL FROM #temp_graph WHERE dest IN (SELECT DISTINCT src FROM #temp_graph);
		-- ignore vertices whose in-degree or out-degree is zero, which cannot appear in any loop.
		DELETE FROM #temp_graph WHERE (src NOT IN (SELECT DISTINCT vertex FROM #temp_link) ) OR (dest NOT IN (SELECT DISTINCT vertex FROM #temp_link) ) ;
		SELECT @var_helper1 = COUNT(*) FROM #temp_graph;
   END
   
   INSERT INTO #temp_graph_backup (src, dest) SELECT src, dest FROM #temp_graph;   
   
   SELECT TOP 1 @var_vertex=src FROM #temp_graph;
   WHILE (@var_vertex IS NOT NULL)
   BEGIN
		SET @var_adj = NULL;
		SET @var_call_stack_index = 0;
		SET @var_vertex_stack_count = 0;
		-- SELECT "Begin vertex", @var_vertex;
		-- SELECT * FROM #temp_graph;
	begin_of_find:
		WHILE (@var_call_stack_index >=0)
		BEGIN
			IF (@var_adj IS NULL)
			BEGIN
				-- first call;
				UPDATE #temp_link SET num = @var_vertex_index_global, in_vertex_stack = 1 WHERE vertex = @var_vertex;
				INSERT #temp_vertex_stack(idx, vertex) VALUES(@var_vertex_stack_count, @var_vertex);
				SET @var_vertex_index = @var_vertex_index_global;
				SET @var_vertex_index_global = @var_vertex_index_global + 1;
				SET @var_vertex_stack_count = @var_vertex_stack_count + 1;
			END
			ELSE
			BEGIN
				-- return from previous call
				UPDATE #temp_link SET num = @var_ret_value WHERE vertex = @var_vertex AND @var_ret_value < num;
				SET @var_adj = NULL;
			END
			
			-- adj is NULL at this point, pick the next available adjacency vertex.
			SELECT TOP 1 @var_adj=dest FROM #temp_graph WHERE src = @var_vertex;
			-- SELECT @var_vertex, @var_vertex_index, @var_adj;
			WHILE (@var_adj IS NOT NULL)
			BEGIN
				-- if there exists one adjacency vertex, visit it.
				-- delete the edge to avoid being revisited.				
				DELETE FROM #temp_graph WHERE src = @var_vertex AND dest = @var_adj;
				SET @var_adj_link_num = NULL;
				SET @var_adj_in_vertex_stack = NULL;
				SELECT @var_adj_link_num = num, @var_adj_in_vertex_stack = in_vertex_stack FROM #temp_link WHERE vertex = @var_adj;
				-- SELECT @var_adj_link_num, @var_adj_in_vertex_stack;
				IF (@var_adj_link_num IS NULL)
				BEGIN
					-- if this adjacency vertex hasn't already visited
					-- save current to call stack
					INSERT #temp_call_stack(idx, vertex, adj, vertex_index) VALUES(@var_call_stack_index, @var_vertex, @var_adj, @var_vertex_index);
					SET @var_call_stack_index = @var_call_stack_index + 1;
					SET @var_vertex = @var_adj;
					SET @var_adj = NULL;
					GOTO begin_of_find;
				END
				ELSE IF (@var_adj_in_vertex_stack IS NOT NULL)
				BEGIN
					UPDATE #temp_link SET num = @var_adj_link_num WHERE vertex = @var_vertex AND @var_adj_link_num < num;
				END
				SET @var_adj = NULL;
				SELECT TOP 1 @var_adj=dest FROM #temp_graph WHERE src = @var_vertex;
				-- SELECT @var_vertex, @var_vertex_index, @var_adj;
			END		
			
							
			-- all the adjacency vertices has already been visited;
				
			-- check if the lowest link of the vertex is the same as that of it's 
			-- vertex_index, if yes a connected component has been found.
			SELECT @var_helper = COUNT(*) FROM #temp_link WHERE vertex = @var_vertex AND num = @var_vertex_index;
			IF (@var_helper > 0)
			BEGIN
				-- SELECT "Find one component", @var_vertex_index;
				-- a strongly connected component has been found.
				SELECT @var_vertex_stack_count = idx FROM #temp_vertex_stack WHERE vertex = @var_vertex;
				INSERT INTO #temp_scc(vertex, component) SELECT vertex, @var_component FROM #temp_vertex_stack WHERE idx >= @var_vertex_stack_count;
				UPDATE #temp_link SET in_vertex_stack = NULL WHERE vertex IN (SELECT DISTINCT vertex FROM #temp_vertex_stack WHERE idx >= @var_vertex_stack_count);
				DELETE FROM #temp_vertex_stack WHERE idx >= @var_vertex_stack_count;
				SET @var_component = @var_component + 1;
			END
			
			-- this call is finished, retun value in ret_value
			-- and recover call stack.
			SELECT @var_ret_value = num FROM #temp_link WHERE vertex = @var_vertex;
			SET @var_call_stack_index = @var_call_stack_index - 1;
			SELECT @var_vertex = vertex, @var_adj = adj, @var_call_stack_index = idx, @var_vertex_index = vertex_index FROM #temp_call_stack WHERE idx = @var_call_stack_index;
			DELETE FROM #temp_call_stack WHERE idx = @var_call_stack_index;
		END
		
		SET @var_vertex = NULL;
		SELECT TOP 1 @var_vertex = src FROM #temp_graph;
   END
   
   -- return the result
   SELECT * FROM #temp_graph_backup, #temp_scc t1, #temp_scc t2 WHERE src = dest OR (src = t1.vertex AND dest = t2.vertex AND t1.component = t2.component);
END 