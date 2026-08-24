# Recommendation Systems Project — Bright Cape Collaboration

BSc team project (4 students), Vrije Universiteit Amsterdam, in collaboration with Bright Cape, a data analytics consultancy (Dec 2024 – March 2025).

Built and compared four recommendation engines (favorite-category, item-based co-occurrence, and user-/item-based collaborative filtering) on a 540,000-observation UK retail transaction dataset, evaluated using precision@5, recall@5, hit-ratio, and NDCG.

**My individual contribution**: the product categorization module — a rule-based classification pipeline (`Bright_Cape_code.Rmd`) that structured ~240,000 retail observations into 34 product categories using a keyword-dictionary approach with priority-weighted matching for products spanning multiple category keywords. This categorized output fed directly into the team's downstream recommendation models.

Data source: UCI Online Retail dataset (not included here).
