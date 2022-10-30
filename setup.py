from pathlib import Path
from typing import List

from setuptools import setup

this_directory = Path(__file__).parent


def get_requirements(filename: str, remove_links=True) -> List[str]:
    """
    lists the requirements to install.
    """
    with open(filename) as f:
        requirements = f.read().splitlines()
    if remove_links:
        for requirement in requirements:
            # git repository url.
            if requirement.startswith("git+"):
                requirements.remove(requirement)
            # subversion repository url.
            if requirement.startswith("svn+"):
                requirements.remove(requirement)
            # mercurial repository url.
            if requirement.startswith("hg+"):
                requirements.remove(requirement)
    return requirements


def get_links(filename: str) -> List[str]:
    """
    gets URL Dependency links.
    """
    links_list = get_requirements(filename, remove_links=False)
    for link in links_list:
        keep_link = False
        already_removed = False
        # git repository url.
        if not link.startswith("git+"):
            if not link.startswith("svn+"):
                if not link.startswith("hg+"):
                    links_list.remove(link)
                    already_removed = True
                else:
                    keep_link = True
                if not keep_link and not already_removed:
                    links_list.remove(link)
                    already_removed = True
            else:
                keep_link = True
            if not keep_link and not already_removed:
                links_list.remove(link)
    return links_list


if __name__ == "__main__":
    use_scm_version = {"root": ".", "relative_to": __file__, "local_scheme": "no-local-version"}
    setup(
        use_scm_version=use_scm_version,
        install_requires=get_requirements("requirements.txt"),
        dependency_links=get_links("requirements.txt")
    )
