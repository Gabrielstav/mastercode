# Call things from here later


# def print_metrics(self):
#     for graph_name, metrics in self.metrics_dict.items():
#         print(f"Metrics for {graph_name}:")
#         for metric_name, metric_value in metrics.items():
#             if callable(metric_value):
#                 if isinstance(metric_value, types.MethodType):
#                     if metric_value.__self__ is None:  # It's a class method
#                         metric_value = metric_value.__func__(type(self))
#                     else:  # It's an instance method
#                         if metric_name == "modularity":
#                             membership = metric_value.__self__.community_multilevel()
#                             metric_value = metric_value.__func__(metric_value.__self__, membership)
#                         else:
#                             metric_value = metric_value.__func__(metric_value.__self__)
#             print(f"  {metric_name}: {metric_value}")
#         print()

# def print_metrics(cls, metrics_data):
#     for graph_name, metrics in metrics_data.items():
#         print(f"Metrics for {graph_name}:")
#         for metric_name, metric_value in metrics.items():
#             print(f"  {metric_name}: {metric_value}")
#         print()

# metrics_data = {}
# for graph_name, graph in graph_dict.items():
#     metrics_data[graph_name] = {}
#     for metric_name, metric_function in metrics.items():
#         if callable(metric_function):
#             metric_value = metric_function(graph)
#         else:
#             metric_value = metric_function
#         metrics_data[graph_name][metric_name] = metric_value
#
# return metrics_data

# def print_metrics(self):
#     for graph_name, metrics in self.metrics_dict.items():
#         print(f"Metrics for {graph_name}:")
#         for metric_name, metric_value in metrics.items():
#             if callable(metric_value):
#                 if isinstance(metric_value, types.MethodType):
#                     if metric_value.__self__ is None:  # It's a class method
#                         metric_value = metric_value.__func__
#                     else:  # It's an instance method
#                         if metric_name == "modularity":
#                             membership = metric_value.__self__.community_multilevel()
#                             metric_value = metric_value.__func__
#                         else:
#                             metric_value = metric_value.__func__
#                         # if metric_name == "modularity":
#                         #     membership = metric_value.__self__.community_multilevel()
#                         #     metric_value = metric_value.__func__(metric_value.__self__, membership)
#                         # else:
#                         #     metric_value = metric_value.__func__(metric_value.__self__)
#                 else:
#                     metric_value = metric_value()
#             print(f"  {metric_name}: {metric_value}")
#         print()

# return cls(graph_dict_or_function, metrics_data)
#
# metrics_data = {}
# for graph_name, graph in graph_dict.items():
#     metrics_data[graph_name] = {}
#     for metric_name, metric_function in metrics.items():
#         metrics_data[graph_name][metric_name] = metric_function(graph)
#
# return cls(graph_dict_or_function, metrics_data)


# for graph_name, graph in graph_dict.items():
#     metrics_data[graph_name] = {}
#     for metric_name, metric_function in metrics.items():
#         try:
#             metric_value = metric_function(graph)
#         except TypeError:
#             # Handle the case when metric_function requires additional arguments
#             if metric_name == "modularity":
#                 membership = graph.community_multilevel().membership
#                 metric_value = metric_function(graph, membership)
#             else:
#                 raise ValueError(f"Unsupported metric '{metric_name}' with additional arguments")
#         metrics_data[graph_name][metric_name] = metric_value

# for graph_name, graph in graph_dict.items():
#     metrics_data[graph_name] = {}
#     for metric_name, metric_function in metrics.items():
#         if callable(metric_function):
#             try:
#                 metric_value = metric_function(graph)
#             except TypeError:
#                 # Handle the case when metric_function requires additional arguments
#                 if metric_name == "modularity":
#                     membership = graph.community_multilevel().membership
#                     metric_value = metric_function(graph, membership)
#                 else:
#                     raise ValueError(f"Unsupported metric '{metric_name}' with additional arguments")
#         else:
#             metric_value = metric_function
#         metrics_data[graph_name][metric_name] = metric_value
