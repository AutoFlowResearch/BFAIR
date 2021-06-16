import React, { useReducer } from 'react';
import { Investigation } from '../../../../../shared/model/investigation.model';
import Card from '../../../../atoms/Card';
import InvestigationFormGeneral from '../InvestigationFormGeneral';
import InvestigationFormStudies from '../InvestigationFormStudies';

const InvestigationFormContent = () => {
  return (
    <div className='investigation-form'>
      <Card>
        <div className='investigation-form__inner'>
          <h3 className='form-title'>Investigation</h3>
          <InvestigationFormGeneral />
          <InvestigationFormStudies/>
        </div>
        <div className='card__footer'>
          <button className='btn--flat'>Cancel</button>
          <button className='btn--primary float--right'>Submit</button>
        </div>
      </Card>
    </div>
  );
};

export default InvestigationFormContent;
