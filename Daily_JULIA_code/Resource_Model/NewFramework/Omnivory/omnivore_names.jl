birmmals_diets = diets[in.(diets[:, 1], Ref(birmmals_names)), :]

plant_items = ["fruit", "other plant parts", "seeds nuts and grains", "algae", "mosses and lichens", "mushrooms"]

# Find species that consume plant items
plant_eaters = unique(birmmals_diets[in.(birmmals_diets[:, 2], Ref(plant_items)), 1])

# Omnivores = plant-eating species that are in predator list but not in herbivore list
omnivore_names = intersect(plant_eaters, predator_names)

carnivore_names = setdiff(predator_names, omnivore_names)
