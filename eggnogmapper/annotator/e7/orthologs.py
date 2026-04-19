"""Ortholog retrieval and donor tracking for v7+ databases.

Provides methods for:
- Finding orthologs from speciation events
- Tracking closest donors for functional annotation
- Filtering orthologs by taxonomic scope
"""

from typing import Dict, List, Optional, Set, Tuple

from .db import EggnogDB
from .tax_scope import LineageFilter
from ..codec import decode_intlist


class OrthologFinder:
    """Find orthologs from speciation events.

    Uses event_index for fast protein -> events lookup,
    then retrieves ortholog IDs from the other side of each event.
    """

    def __init__(
        self,
        db: EggnogDB,
        lineage_filter: Optional[LineageFilter] = None,
    ):
        """Initialize ortholog finder.

        Args:
            db: EggnogDB instance
            lineage_filter: Optional LineageFilter for tax_scope filtering
        """
        self.db = db
        self.lineage_filter = lineage_filter

    def find_orthologs(
        self,
        protein_id: int,
        seed_taxid: Optional[str] = None,
    ) -> List[int]:
        """Find all orthologs for a protein.

        Args:
            protein_id: Integer protein ID
            seed_taxid: Seed species taxid (for auto tax_scope)

        Returns:
            List of ortholog protein IDs (deduplicated)
        """
        # Get event IDs from index
        event_ids = self.db.get_events_for_protein(protein_id)
        if not event_ids:
            return []

        # Fetch event data
        events = self.db.get_events_bulk(event_ids)

        # Collect orthologs from opposite side of each event
        orthologs = set()
        for event_id, event in events.items():
            if protein_id in event["side1"]:
                ortholog_ids = event["side2"]
            elif protein_id in event["side2"]:
                ortholog_ids = event["side1"]
            else:
                continue

            orthologs.update(ortholog_ids)

        ortholog_list = list(orthologs)

        # Apply taxonomic scope filter if configured
        if self.lineage_filter:
            ortholog_list = self.lineage_filter.filter_ortholog_ids(
                ortholog_list,
                seed_taxid=seed_taxid,
            )

        return ortholog_list

    def find_orthologs_with_events(
        self,
        protein_id: int,
        seed_taxid: Optional[str] = None,
    ) -> List[Tuple[int, dict]]:
        """Find orthologs with their event metadata.

        Args:
            protein_id: Integer protein ID
            seed_taxid: Seed species taxid (for auto tax_scope)

        Returns:
            List of (ortholog_id, event_dict) tuples
        """
        event_ids = self.db.get_events_for_protein(protein_id)
        if not event_ids:
            return []

        events = self.db.get_events_bulk(event_ids)

        results = []
        for event_id, event in events.items():
            if protein_id in event["side1"]:
                ortholog_ids = event["side2"]
            elif protein_id in event["side2"]:
                ortholog_ids = event["side1"]
            else:
                continue

            # Apply tax scope filter
            if self.lineage_filter:
                ortholog_ids = self.lineage_filter.filter_ortholog_ids(
                    ortholog_ids,
                    seed_taxid=seed_taxid,
                )

            for oid in ortholog_ids:
                results.append((oid, event))

        return results


class DonorTracker:
    """Track closest donors for functional annotation.

    Uses event metadata (og_lca, ev_lca) to identify the closest
    taxonomic donor for each annotation type.
    """

    def __init__(self, db: EggnogDB):
        """Initialize donor tracker.

        Args:
            db: EggnogDB instance
        """
        self.db = db

    def get_closest_donor(
        self,
        protein_id: int,
        annotation_type: str,
    ) -> Optional[dict]:
        """Find the closest donor for a specific annotation type.

        Args:
            protein_id: Integer protein ID
            annotation_type: Annotation field (e.g., 'kegg_ko', 'gos')

        Returns:
            Dict with donor info, or None if no donor found
        """
        # Get events for protein
        event_ids = self.db.get_events_for_protein(protein_id)
        if not event_ids:
            return None

        events = self.db.get_events_bulk(event_ids)

        # Find ortholog IDs
        all_orthologs = set()
        event_for_ortholog = {}

        for event_id, event in events.items():
            if protein_id in event["side1"]:
                ortholog_ids = event["side2"]
            elif protein_id in event["side2"]:
                ortholog_ids = event["side1"]
            else:
                continue

            for oid in ortholog_ids:
                all_orthologs.add(oid)
                if oid not in event_for_ortholog:
                    event_for_ortholog[oid] = event

        if not all_orthologs:
            return None

        # Fetch annotations for orthologs
        annotations = self.db.get_protein_annotations_bulk(list(all_orthologs))

        # Find first ortholog with the annotation
        for oid, annot in annotations.items():
            if annot and annot.get(annotation_type):
                event = event_for_ortholog.get(oid)
                return {
                    "donor_id": oid,
                    "value": annot[annotation_type],
                    "cluster": event.get("name") if event else None,
                    "og": event.get("og") if event else None,
                    "ev_lca": event.get("ev_lca") if event else None,
                }

        return None

    def get_all_donors(
        self,
        protein_id: int,
    ) -> Dict[str, dict]:
        """Get closest donors for all annotation types.

        Args:
            protein_id: Integer protein ID

        Returns:
            Dict mapping annotation_type -> donor info
        """
        # Annotation fields to track
        annotation_fields = [
            "pname", "gos", "kegg_ko", "kegg_ec", "kegg_pathway",
            "kegg_module", "pfam", "kegg_cog", "kegg_cazy",
        ]

        donors = {}
        for field in annotation_fields:
            donor = self.get_closest_donor(protein_id, field)
            if donor:
                donors[field] = donor

        return donors
