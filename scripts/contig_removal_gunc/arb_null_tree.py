from gtdblib.util.bio.accession import canonical_gid

from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged


def main():
    df_css = AggregateMaxCssLevelMerged().output().maybe_read_cached()
    gids = {canonical_gid(x) for x in df_css.index}

    print(len(gids), list(gids)[:10])

    with open('/tmp/null_tree.tree', 'w') as f:
        f.write(f'({",".join(sorted(gids))});')


    return

if __name__ == '__main__':
    main()