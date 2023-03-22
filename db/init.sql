\connect surechembl;

create schema if not exists sc;
set search_path to sc;

create extension if not exists rdkit;

drop table if exists compound cascade;
drop table if exists tmp_compound cascade;

drop table if exists patent cascade;
drop table if exists tmp_patent cascade;

drop table if exists field cascade;
drop table if exists tmp_field cascade;

drop table if exists field_freq cascade;
drop table if exists tmp_field_freq cascade;

create table compound (
    id int,
    smiles mol not null,
    mfp bfp not null,
    primary key(id)
);

create index smi_idx on compound using gist(smiles);
create index mfp_idx on compound using gist(mfp);

create table patent (
    id serial,
    num varchar(32) unique not null,
    pub_date timestamp not null,
    primary key(id)
);

create table field (
    id int,
    details varchar(64) not null,
    primary key(id)
);

insert into field values (1, 'Description'), 
    (2, 'Claims'), (3, 'Abstract'), (4, 'Title'),
    (5, 'Image (for patents after 2007)'),
    (6, 'MOL Attachment (US patents after 2007)');

create table field_freq (
    compound_id int,
    patent_id int,
    field_id int,
    freq int not null,
    foreign key(compound_id) references compound(id),
    foreign key(patent_id) references patent(id),
    foreign key(field_id) references field(id),
    unique(compound_id, patent_id, field_id)
);
