"""empty message

Revision ID: 0542afe45e5c
Revises: 
Create Date: 2021-02-24 10:07:35.211744

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '0542afe45e5c'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_constraint('Methodology_Sample_id_fkey', 'Methodology', type_='foreignkey')
    op.create_foreign_key(None, 'Methodology', 'Sample', ['Sample_id'], ['id'], ondelete='SET NULL')
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_constraint(None, 'Methodology', type_='foreignkey')
    op.create_foreign_key('Methodology_Sample_id_fkey', 'Methodology', 'Sample', ['Sample_id'], ['id'])
    # ### end Alembic commands ###
